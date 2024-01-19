import f90nml
import os
import numpy as np


class ramses_sim:
    class param_list(dict):
        # pass

        def print_dict(self):
            for key, value in self.__dict__.items():
                print(key, value)

        def __getattr__(self, item):
            return super().__getitem__(item)

        def __setattr__(self, item, value):
            return super().__setitem__(item, value)

    def __init__(self, path, nml=None):
        self.path = path

        self.namelist = self.param_list()
        nml_params = get_nml_params(path, name=nml)
        for key, value in nml_params.items():
            setattr(self.namelist, key, value)

        snaps, snap_numbers = get_snaps(path)
        self.snaps = snaps
        self.snap_numbers = snap_numbers

        first_snap = snap_numbers[0]
        first_info_params = get_info_params(path, first_snap)

        last_info_params = get_info_params(path, snap_numbers[-1])

        self.aexp_stt = np.float64(first_info_params["aexp"])
        self.aexp_end = np.float64(last_info_params["aexp"])

        self.ncpu = int(first_info_params["ncpu"])
        self.ndim = int(first_info_params["ndim"])
        self.levelmin = int(first_info_params["levelmin"])
        self.levelmax = int(first_info_params["levelmax"])

        self.cosmo = self.param_list()
        self.cosmo.boxlen = np.float64(first_info_params["boxlen"])
        self.cosmo.H0 = np.float64(first_info_params["H0"])
        self.cosmo.Omega_m = np.float64(first_info_params["omega_m"])
        self.cosmo.Omega_l = np.float64(first_info_params["omega_l"])
        self.cosmo.Omega_b = np.float64(first_info_params["omega_b"])
        self.cosmo.Omega_k = np.float64(first_info_params["omega_k"])
        self.cosmo.unit_d = np.float64(first_info_params["unit_d"])
        self.cosmo.unit_l = np.float64(first_info_params["unit_l"])
        self.cosmo.unit_t = np.float64(first_info_params["unit_t"])

        hydroparams = read_hydro_file(path, first_snap)
        if hydroparams != None:
            self.hydro = self.param_list()
            for key, value in hydroparams.items():
                setattr(self.hydro, key, value)

    def get_snap_exps(self, snap_nbs=None):
        if not "aexps" in os.listdir(self.path):
            if snap_nbs is None:
                snap_nbs = self.snap_numbers

            aexps = np.zeros(len(snap_nbs), dtype="f4")

            for isnap, snap_nb in enumerate(snap_nbs):
                aexps[isnap] = get_info_params(self.path, snap_nb)["aexp"]

            np.save(os.path.join(self.path, "aexps"), aexps)

        else:
            aexps = np.load(os.path.join(self.path, "aexps"))

        self.aexps = aexps

        return aexps


def get_snaps(path, SIXDIGITS=False):
    snaps = np.sort([x for x in os.listdir(os.path.join(path)) if "output_" in x])
    snap_numbers = np.array([int(x.split("_")[-1]) for x in snaps])

    info_present = np.zeros(len(snap_numbers), dtype=bool)

    for isnap, snap_nb in enumerate(snap_numbers):
        try:
            get_info_params(path, snap_nb, SIXDIGITS=SIXDIGITS)
            info_present[isnap] = True
        except:
            FileNotFoundError

    snaps = snaps[info_present]
    snap_numbers = snap_numbers[info_present]

    arg = np.argsort(snap_numbers)

    return (snaps[arg], snap_numbers[arg])


def get_nml_params(path, name=None):
    if name is None:
        nml = [f for f in os.listdir(path) if ".nml" in f and "~" not in f]
        # print(nml)
        assert len(nml) == 1, "More than one namelist file found"
        nml = nml[0]

    else:
        nml = name

    namelist = f90nml.read(os.path.join(path, nml))

    for cat in namelist.keys():
        for key in namelist[cat].keys():
            # print(cat, key, namelist[cat][key])
            if (
                type(namelist[cat][key]) == list
                and type(namelist[cat][key][0]) == str
                and "!" in namelist[cat][key][0]
            ):
                # print("caught")
                namelist[cat][key] = namelist[cat][key][0].split("!")[0].strip()

    return namelist


def read_info_file(fname):
    infos = {}

    with open(fname, "r") as f:
        for il, line in enumerate(f):
            # if line[0] != '#' or line[0] != ' ':
            if "=" in line:
                # print('|',line,'|')
                key, value = line.split("=")
                key = key.strip()
                value = value.strip()
                if not "." in value:
                    infos[key] = int(value)
                else:
                    infos[key] = float(value)

                infos[key] = value

            if il == 17:
                break
    return infos


def get_info_params(path, snap, SIXDIGITS=False):
    if not SIXDIGITS:
        snap_str = f"{snap:05d}"
    else:
        snap_str = f"{snap:06d}"

    infos = {}

    with open(
        os.path.join(path, f"output_{snap_str}", f"info_{snap_str}.txt"), "r"
    ) as f:
        for il, line in enumerate(f):
            # if line[0] != '#' or line[0] != ' ':
            if "=" in line:
                # print('|',line,'|')
                key, value = line.split("=")
                key = key.strip()
                value = value.strip()
                if not "." in value:
                    infos[key] = int(value)
                else:
                    infos[key] = float(value)

                infos[key] = value

            if il == 17:
                break

    return infos


def read_hydro_file(path, snap, SIXDIGITS=False):
    if not SIXDIGITS:
        fname = os.path.join(path, f"output_{snap:05d}", "hydro_file_descriptor.txt")
    else:
        fname = os.path.join(path, f"output_{snap:06d}", "hydro_file_descriptor.txt")
    with open(fname, "r") as src:
        hydro = {}

        for line in src:
            # print(line)
            if line[0] == "#":
                continue
            elif "=" in line:
                key, value = line.split("=")
                key = key.strip()
                value = value.strip()

                hydro[key] = int(value)
            elif ":" in line:
                key, value = line.split(":")
                key = key.strip()
                value = value.strip()
            elif "," in line:
                key, value, dtype = line.split(",")
                key = key.strip()
                value = value.strip()

                hydro[key] = value

    return hydro
