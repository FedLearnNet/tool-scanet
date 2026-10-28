from pyfedappwrap.engine.runtime import FedDBEngine

from app import MyAPP

engine = FedDBEngine()
engine.register(MyAPP())

if __name__ == "__main__":
    engine.start()
    engine.wait_until_stop()
