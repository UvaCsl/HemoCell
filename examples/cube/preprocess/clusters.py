import abc


class Cluster(abc.ABC):
    def __init__(self, cpu_per_node):
        self.cpu_per_node = cpu_per_node

    @abc.abstractmethod
    def partition(self, node_count):
        ...


class Lisa(Cluster):
    def __init__(self):
        super().__init__(16)
        self.environment_script = 'lisa_env.sh'
        self.energy_collection = ''

    def partition(self, node_count):
        return 'thin'


class Snellius(Cluster):
    def __init__(self):
        super().__init__(64)
        self.environment_script = 'snellius_env.sh'
        self.energy_collection = ''

    def partition(self, node_count):
        return 'thin'


class SuperMUC(Cluster):
    def __init__(self):
        super().__init__(48)
        self.environment_script = 'lrz_env.sh'
        self.energy_collection = ''

    def partition(self, node_count):
        if node_count < 16:
            return 'micro'
        if node_count < 768:
            return 'general'
        return 'large'


class Archer2(Cluster):
    def __init__(self):
        super().__init__(128)
        self.environment_script = 'archer2_env.sh'
        self.energy_collection = 'sacct -j ${SLURM_JOB_ID} --format=ALL' # To convert to [kWh] we can multiply the energy in  [J] joules by 2.78e-7

    # Archer2 requires a "quality of service" flag in addition to the partition
    # specification to match with the requested job specifications.
    def qos(self, node_count):
        if node_count <= 1024:
            return 'standard'
        if node_count <= 5860:
            return 'largescale'

    def partition(self, node_count):
        return 'standard'


def cluster_from_string(string):
    if string == 'snellius':
        return Snellius()
    if string == 'supermuc':
        return SuperMUC()
    if string == 'lisa':
        return Lisa()
    if string == 'archer2':
        return Archer2()
