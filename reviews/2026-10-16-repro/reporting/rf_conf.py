import importlib.util, sys
spec = importlib.util.spec_from_file_location("cc", "./src/robin/classification_config.py")
cc = importlib.util.module_from_spec(spec); spec.loader.exec_module(cc)
for name in ["Sturgeon","NanoDX","PanNanoDX","Random Forest"]:
    key = name.lower()
    print(name, repr(key), key in cc.CLASSIFIER_CONFIDENCE_THRESHOLDS)
for v in [0.70, 0.82]:
    print("report RF", v, cc.get_confidence_status("random forest", v), "| GUI RF", cc.get_confidence_ui_tier("random_forest", v*100), cc.get_confidence_level("random_forest", v*100))
