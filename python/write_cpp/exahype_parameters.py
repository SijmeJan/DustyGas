class ExaHyPE_parameters():
    def __init__(self, filename):
        # Read exahype file
        f = open(filename, "r")
        lines = f.readlines()
        f.close()

        # Remove comments between /* and */
        start_comment = []
        end_comment = []
        for i in range(0, len(lines)):
            if (lines[i].find('/*') != -1):
                start_comment.append(i)
            if (lines[i].find('*/') != -1):
                end_comment.append(i)
        for n in range(len(start_comment)-1, -1, -1):
            for i in range(0, end_comment[n] - start_comment[n] + 1):
                lines.pop(start_comment[n])

        # Parameters to be read from exahype file
        self.n_vars = 0
        self.patch_size = 0
        self.order = None
        self.output_dir = None
        self.project_name = None
        self.solver_name = None
        self.solver_type = None
        self.boundary_name = None
        self.offset_x = None
        self.offset_y = None
        self.size_x = None
        self.size_y = None
        self.cell_size = None

        found_nvar = False
        for line in lines:
            if (line.find('exahype-project ') != -1):
                self.project_name = line.lstrip().split()[-1]
            if (line.find('solver ') != -1):
                self.solver_name = line.lstrip().split()[-1]
                self.solver_type = line.lstrip().split()[-2]
            if (line.find('variables const') != -1 and found_nvar == False):
                self.n_vars = int(line.lstrip().split()[-1])
                found_nvar = True
            if (line.find('output-directory') != -1):
                self.output_dir = line.lstrip().split()[-1]
            if (line.find('width') != -1):
                self.size_x = float(line.lstrip().split()[-2][:-1])
                self.size_y = float(line.lstrip().split()[-1])
            if (line.find('offset') != -1):
                self.offset_x = float(line.lstrip().split()[-2][:-1])
                self.offset_y = float(line.lstrip().split()[-1])
            if (line.find('patch-size const') != -1):
                self.patch_size = int(line.lstrip().split()[-1])
            if (line.find('order const') != -1):
                self.order = int(line.lstrip().split()[-1])
            if (line.find('maximum-mesh-size') != -1):
                self.cell_size = float(line.lstrip().split()[-1])
