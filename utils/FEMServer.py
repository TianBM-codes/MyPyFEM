from main import MyPyFEM, AnalyseType
import pathlib

class FEMServer:
    _instance = None

    def __init__(self, input_file):
        self.fem = MyPyFEM(
            pathlib.Path(input_file),
            check_model=False,
            plot_stiff=False,
            AnaType=AnalyseType.AsServer
        )

    @classmethod
    def initialize(cls, input_file):
        if cls._instance is None:
            cls._instance = cls(input_file)
        return cls._instance

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            raise RuntimeError("FEMServer not initialized! Call FEMServer.initialize() first.")
        return cls._instance