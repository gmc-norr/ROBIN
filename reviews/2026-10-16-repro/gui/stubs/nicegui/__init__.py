class _Any:
    def __getattr__(self, name):
        return _Any()
    def __call__(self, *a, **k):
        if len(a) == 1 and callable(a[0]) and not k:
            return a[0]
        return _Any()
ui = _Any()
app = _Any()
run = _Any()
class _Ctx: pass
context = _Any()
