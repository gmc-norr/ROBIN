TIMERS = []
DISCONNECT = []
class El:
    def __init__(self, *a, **k):
        self.text = a[0] if a and isinstance(a[0], str) else ""
        self.options = a[0] if a and isinstance(a[0], dict) else {}
        self.value = k.get("value")
    def __getattr__(self, name):
        def f(*a, **k):
            return self
        return f
    def set_text(self, t): self.text = t
    def __enter__(self): return self
    def __exit__(self, *a): return False
    def update(self): pass
class _Client:
    def on_disconnect(self, cb): DISCONNECT.append(cb)
    def on_connect(self, cb): pass
class _Context:
    client = _Client()
class _UI:
    context = _Context()
    def timer(self, interval, cb, once=False, active=True, immediate=True):
        t = El(); t.interval = interval; t.cb = cb; t.once = once; t.active = active
        def deactivate(): t.active = False
        t.deactivate = deactivate
        TIMERS.append(t); return t
    def refreshable(self, f): return f
    def __getattr__(self, name):
        return lambda *a, **k: El(*a, **k)
ui = _UI()
class _Storage:
    user = {}
    general = {}
class _App:
    storage = _Storage()
app = _App()
run = El()
