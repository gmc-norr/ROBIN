#!/bin/bash
# Reproduce: update_clinvar_if_newer replaces clinvar.vcf.gz but keeps the old .tbi
set -u
D=$(mktemp -d)
mkdir -p "$D"
cd "$D" || exit 1
python3 - <<'EOF'
import random
hdr = "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"x\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
random.seed(1)
def write(fn, n, offset):
    with open(fn, "w") as f:
        f.write(hdr)
        for i in range(n):
            pos = 1000 + i * 50 + offset
            f.write(f"1\t{pos}\tv{i}\tA\tG\t.\t.\tCLNSIG=Pathogenic_{'x'*random.randint(0,200)}\n")
write("old.vcf", 200000, 0)
write("new.vcf", 260000, 7)
EOF
bcftools view -Oz -o clinvar.vcf.gz old.vcf && bcftools index -f -t clinvar.vcf.gz
echo "== old file, old index:"; bcftools view -H -r 1:5000000-5000200 clinvar.vcf.gz | cut -f1-3
bcftools view -Oz -o new.vcf.gz new.vcf && mv -f new.vcf.gz clinvar.vcf.gz
touch -d '1 hour ago' clinvar.vcf.gz.tbi
echo "== new file, STALE old index:"; bcftools view -H -r 1:5000000-5000200 clinvar.vcf.gz 2>&1 | cut -f1-3 | head
echo "   exit=$?"
cp -f clinvar.vcf.gz fresh.vcf.gz && bcftools index -f -t fresh.vcf.gz
echo "== new file, fresh index (expected):"; bcftools view -H -r 1:5000000-5000200 fresh.vcf.gz | cut -f1-3
