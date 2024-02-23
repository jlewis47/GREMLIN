import f90nml
import os
import numpy as np
from astropy.cosmology import FlatLambdaCDM


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

    def __init__(self, path, nml=None, output_path=None, info_path=None):
        self.path = path
        self.rel_output_path = output_path
        self.rel_info_path = info_path

        if output_path is None:
            self.output_path = self.path
        else:
            self.output_path = os.path.join(self.path, self.rel_output_path)

        if info_path is None:
            self.info_outputs = True
            self.info_path = self.path
        else:
            self.info_outputs = False
            self.info_path = os.path.join(self.path, info_path)

        self.namelist = self.param_list()
        nml_params = get_nml_params(path, name=nml)
        for key, value in nml_params.items():
            setattr(self.namelist, key, value)

        snaps, snap_numbers = self.get_snaps(path)
        self.snaps = snaps
        self.snap_numbers = snap_numbers

        first_snap = snap_numbers[0]
        last_snap = snap_numbers[-1]
        first_info_params = get_info_params(
            self.info_path, first_snap, outputs=self.info_outputs
        )

        last_info_params = get_info_params(
            self.info_path, last_snap, outputs=self.info_outputs
        )

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

        try:
            hydroparams = read_hydro_file(path, last_snap)
            if hydroparams != None:
                self.hydro = self.param_list()
                for key, value in hydroparams.items():
                    setattr(self.hydro, key, value)
        except FileNotFoundError:
            print("No hydro file found")

    def init_cosmo(self):
        self.cosmo_model = FlatLambdaCDM(
            H0=self.cosmo.H0,
            Om0=self.cosmo.Omega_m,
            Ob0=self.cosmo.Omega_b,
            Tcmb0=2.725,
        )

    def get_snap_exps(self, snap_nbs=None, param_save=True):
        save = False

        if not "aexps" in os.listdir(self.path):
            if snap_nbs is None:
                snap_nbs = self.snap_numbers
                save = True
            elif type(snap_nbs) == int:
                snap_nbs = [snap_nbs]

            aexps = np.zeros(len(snap_nbs), dtype="f4")

            for isnap, snap_nb in enumerate(snap_nbs):
                aexps[isnap] = get_info_params(
                    self.info_path, snap_nb, outputs=self.info_outputs
                )["aexp"]

            if param_save and save:
                np.save(os.path.join(self.path, "aexps"), aexps)

        else:
            aexps = np.load(os.path.join(self.path, "aexps"))

        self.aexps = aexps

        return aexps

    def get_snap_times(self, snap_nbs=None, param_save=True):

        save = False

        if not hasattr(self, "cosmo_model"):
            self.init_cosmo()

        if snap_nbs is None:
            snap_nbs = self.snap_numbers
            save = True
        elif type(snap_nbs) == int:
            snap_nbs = [snap_nbs]

        if not hasattr(self, "aexps"):
            self.get_snap_exps(snap_nbs=snap_nbs, param_save=param_save)

        if not "times" in os.listdir(self.path):

            times = self.cosmo_model.age(1.0 / self.aexps - 1.0).value * 1e3  # Myr

            if param_save and save:
                np.save(os.path.join(self.path, "times"), times)

        else:

            times = np.load(os.path.join(self.path, "times"))

        self.times = times

        return times

    def get_snaps(self, path, SIXDIGITS=False):
        snaps = np.sort([x for x in os.listdir(path) if "output_" in x])
        snap_numbers = np.array([int(x.split("_")[-1]) for x in snaps])

        info_present = np.zeros(len(snap_numbers), dtype=bool)

        for isnap, snap_nb in enumerate(snap_numbers):
            try:
                get_info_params(
                    self.info_path,
                    snap_nb,
                    SIXDIGITS=SIXDIGITS,
                    outputs=self.info_outputs,
                )
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


def get_info_params(path, snap, SIXDIGITS=False, outputs=True):
    if not SIXDIGITS:
        snap_str = f"{snap:05d}"
    else:
        snap_str = f"{snap:06d}"

    infos = {}

    if outputs:
        info_fpath = os.path.join(path, f"output_{snap_str}", f"info_{snap_str}.txt")
    else:
        info_fpath = os.path.join(path, f"info_{snap_str}.txt")

    with open(info_fpath, "r") as f:
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
                hydro[key] = value
            elif "," in line:
                key, value, dtype = line.split(",")
                key = key.strip()
                value = value.strip()

                hydro[key] = value

    return hydro
