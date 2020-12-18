import os
from time import time
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import pyemu

par_tags = {"total_erupted_mass", "column_height", "ellipse_major_axis", "ellipse_minor_axis", "diffusion_coef",
            "wind_speed", "wind_direction", "tgsd_mean", "tgsd_sigma"}


def write_config_tpl_file(config_filename):
    """write a PEST-style template file for the config file. This
	template file is used by PEST(++) tools to write the desired parameter
	values into a file that the simulation can use.

	Harvest the values in the config value and return them in a dataframe

	"""
    lines = open(config_filename, 'r').readlines()
    tpl_filename = config_filename + ".tpl"
    par_names, par_vals = [], []
    with open(tpl_filename, 'w') as f:
        f.write("ptf ~\n")
        for line in lines:
            if line.lower().split("=")[0].strip() in par_tags:
                par_name = line.lower().split("=")[0].strip()
                par_names.append(par_name)
                try:
                    par_val = float(line.lower().split("=")[1].split()[0])
                    par_vals.append(par_val)
                except Exception as e:
                    raise Exception("error casting {0} for float".format(line.lower().split("=")[1].split()[0]))

                line = line.split("=")[0] + "=  ~  " + par_name + "  ~ \n"
            f.write(line)
    df = pd.DataFrame({"parnme": par_names, "parval1": par_vals}, index=par_names)
    return tpl_filename, df


def write_output_ins_file(output_filename):
    """write a PEST-style instruction file for the output file.
	This will instruct PEST(++) tools on how to read the outputs from the
	simulation that correspond to the observations

	Harvest the output values in the output file and return them in a
	dataframe

	"""
    lines = open(output_filename, 'r').readlines()
    ins_filename = output_filename + ".ins"
    obs_names, obs_vals = [], []
    with open(ins_filename, 'w') as f:
        f.write("pif ~\n")
        f.write("l1 \n")  # skip the header line
        for line in lines[1:]:  # skip the header line
            raw = line.strip().split()
            obs_name = "obs_{0}_x:{1}_y:{2}".format(raw[0], raw[1], raw[2])
            assert obs_name not in obs_names, obs_name
            obs_names.append(obs_name)
            try:
                obs_val = float(raw[-1])
                obs_vals.append(obs_val)
            except Exception as e:
                raise Exception("error casting {0} to float".format(raw[-1]))
            f.write("l1 w w w !{0}! \n".format(obs_name))

    df = pd.DataFrame({"obsnme": obs_names, "obsval": obs_vals}, index=obs_names)
    return ins_filename, df


def setup_pulu():
    # run the model one
    shutil.copy2(os.path.join("pulu", "org_config_pulu.py"), os.path.join("pulu", "config_pulu.py"))
    pyemu.os_utils.run("python main.py", cwd="pulu")

    # prepare the interface files
    output_filename = os.path.join("pulu", "pulu_sim.dat")
    ins_filename, obs_df = write_output_ins_file(output_filename)
    config_filename = os.path.join(os.path.join("pulu", "config_pulu.py"))
    tpl_filename, par_df = write_config_tpl_file(config_filename)

    # create a control file instance
    pst = pyemu.Pst.from_io_files(tpl_files=tpl_filename, in_files=config_filename, ins_files=ins_filename,
                                  out_files=output_filename, pst_path=".")

    # set the observation values
    pst.observation_data.loc[obs_df.obsnme, "obsval"] = obs_df.obsval

    # set the parameter information for the pulu pars
    par = pst.parameter_data
    # start by fixing all the parameters in the control file
    par.loc[:, "partrans"] = "fixed"

    par.loc["total_erupted_mass", "parlbnd"] = 1e10
    par.loc["total_erupted_mass", "parubnd"] = 1e12
    par.loc["total_erupted_mass", "parval1"] = 1e11
    par.loc["total_erupted_mass", "partrans"] = "log"

    par.loc["column_height", "parlbnd"] = 8000
    par.loc["column_height", "parubnd"] = 20000
    par.loc["column_height", "parval1"] = 15000
    par.loc["column_height", "partrans"] = "log"

    par.loc["diffusion_coef", "parlbnd"] = 1000
    par.loc["diffusion_coef", "parubnd"] = 10000
    par.loc["diffusion_coef", "parval1"] = 3000
    par.loc["diffusion_coef", "partrans"] = "log"

    # it looks like mean grainsize is already in log space?
    par.loc["tgsd_mean", "parlbnd"] = -5
    par.loc["tgsd_mean", "parubnd"] = +5
    par.loc["tgsd_mean", "parval1"] = -1.99
    par.loc["tgsd_mean", "partrans"] = "none"
    par.loc["tgsd_mean", "parchglim"] = "relative"

    par.loc["ellipse_major_axis", "parlbnd"] = 3000 
    par.loc["ellipse_major_axis", "parubnd"] = 12000
    par.loc["ellipse_major_axis", "parval1"] = 7000
    par.loc["ellipse_major_axis", "partrans"] = "log"

    par.loc["ellipse_minor_axis", "parlbnd"] = 1000
    par.loc["ellipse_minor_axis", "parubnd"] = 7000
    par.loc["ellipse_minor_axis", "parval1"] = 3000
    par.loc["ellipse_minor_axis", "partrans"] = "log"


    # par.loc["ellipse_major_axis", "parlbnd"] = 5000
    # par.loc["ellipse_major_axis", "parubnd"] = 15000
    # par.loc["ellipse_major_axis", "parval1"] = 8000
    # par.loc["ellipse_minor_axis", "parval1"] = 5000

    # # tie the minor axis par to the major so they function
    # # as a single radius
    # par.loc["ellipse_major_axis", "partrans"] = "none"
    # #par.loc["ellipse_minor_axis", "partrans"] = "tied"
    # par.loc["ellipse_minor_axis", "partied"] = "ellipse_major_axis"
  
    par.loc["wind_speed", "parlbnd"] = 0.1
    par.loc["wind_speed", "parubnd"] = 5
    par.loc["wind_speed", "parval1"] = 1
    par.loc["wind_speed", "partrans"] = "log"

    par.loc["wind_direction", "parlbnd"] = 195
    par.loc["wind_direction", "parubnd"] = 225
    par.loc["wind_direction", "parval1"] = 210
    par.loc["wind_direction", "partrans"] = "log"

    # set the observed thickness values
    pst.try_parse_name_metadata()  # this sets individual columns for northing and easting
    obs = pst.observation_data
    # read the observed data
    df = pd.read_csv(os.path.join("pulu", "pulu_grainsize.csv"))
    df.columns = df.columns.map(lambda x: x.lower().replace(" ", "_"))

    # for each observation...
    for oname, oe, on in zip(obs.obsnme, obs.e, obs.n):
        # find the location based on easting and northing (not the relative, the absolute)
        v = df.loc[df.apply(lambda x: x.easting == float(oe) and x.northing == float(on), axis=1), "thickness"].values[
            0]
        # set the observed value
        obs.loc[oname, "obsval"] = v

    # set the model forward run command
    pst.model_command = "python main.py"

    # run pestpp-glm for a single one-off forward run
    pst.control_data.noptmax = 0
    pst.write(os.path.join("pulu", "pulu.pst"))
    pyemu.os_utils.run("pestpp-glm pulu.pst", cwd="pulu")


def pulu_run_prior_mc():
    """run a prior-based monte carlo in parallel locally"""
    pst = pyemu.Pst(os.path.join("pulu", "pulu.pst"))
    pst.control_data.noptmax = -1
    pst.pestpp_options["ies_num_reals"] = 1000

    pst.write(os.path.join("pulu", "pulu_run.pst"))
    # pyemu.os_utils.run("pestpp-ies pulu_run.pst")
    pyemu.os_utils.start_workers("pulu", "pestpp-ies", "pulu_run.pst", num_workers=10,
                                 master_dir="pulu_prior_mc_master")


def pulu_run_ies():
    """run pestpp-ies in parallel locally"""
    pst = pyemu.Pst(os.path.join("pulu", "pulu.pst"))
    pst.control_data.noptmax = 5
    pst.pestpp_options["ies_num_reals"] = 50

    pst.write(os.path.join("pulu", "pulu_run.pst"))
    # pyemu.os_utils.run("pestpp-ies pulu_run.pst")
    pyemu.os_utils.start_workers("pulu", "pestpp-ies", "pulu_run.pst", num_workers=10,
                                 master_dir="pulu_ies_master")


def pulu_run_glm():
    """run pestpp-glm in parallel locally"""
    pst = pyemu.Pst(os.path.join("pulu", "pulu.pst"))
    pst.control_data.noptmax = 5
    pst.pestpp_options["glm_num_reals"] = 100

    pst.write(os.path.join("pulu", "pulu_run.pst"))
    # pyemu.os_utils.run("pestpp-ies pulu_run.pst")
    pyemu.os_utils.start_workers("pulu", "pestpp-glm", "pulu_run.pst", num_workers=10,
                                 master_dir="pulu_glm_master")


def pulu_plot_ies_results():
    """plot the pestpp-ies results"""
    m_d = "pulu_ies_master"
    pst = pyemu.Pst(os.path.join(m_d, "pulu_run.pst"))
    pr_pe = pd.read_csv(os.path.join(m_d, "pulu_run.0.par.csv"), index_col=0)
    pr_oe = pd.read_csv(os.path.join(m_d, "pulu_run.0.obs.csv"), index_col=0)
    pt_pe = pd.read_csv(os.path.join(m_d, "pulu_run.{0}.par.csv".
                                     format(pst.control_data.noptmax)), index_col=0)
    pt_oe = pd.read_csv(os.path.join(m_d, "pulu_run.{0}.obs.csv".
                                     format(pst.control_data.noptmax)), index_col=0)

    pr_pe = pr_pe.loc[:, pst.adj_par_names]
    pt_pe = pt_pe.loc[:, pst.adj_par_names]

    pyemu.plot_utils.ensemble_helper({"0.5": pr_pe, "b": pt_pe}, bins=20,
                                     filename=os.path.join(m_d, "pulu_ies_par_summary.pdf"))
    # plt.show()
    obs = pst.observation_data
    pyemu.plot_utils.ensemble_helper({"0.5": pr_oe, "b": pt_oe},
                                     deter_vals=obs.obsval.to_dict(),
                                     bins=20,
                                     filename=os.path.join(m_d, "pulu_ies_obs_summary.pdf"))
    pyemu.plot_utils.ensemble_res_1to1(pst=pst, ensemble={"0.5": pr_oe, "b": pt_oe},
                                       filename=os.path.join(m_d, "pulu_ies_obs_vs_sim.pdf"))
    plt.show()


def pulu_plot_glm_results():
    """plot the pestpp-glm results"""
    m_d = "pulu_glm_master"
    pst = pyemu.Pst(os.path.join(m_d, "pulu_run.pst"))

    pt_pe = pd.read_csv(os.path.join(m_d, "pulu_run.post.paren.csv"), index_col=0)
    pt_oe = pd.read_csv(os.path.join(m_d, "pulu_run.post.obsen.csv"), index_col=0)

    pt_pe = pt_pe.loc[:, pst.adj_par_names]

    pyemu.plot_utils.ensemble_helper({"b": pt_pe}, bins=40,
                                     filename=os.path.join(m_d, "pulu_glm_par_summary.pdf"))
    # plt.show()
    obs = pst.observation_data
    pyemu.plot_utils.ensemble_helper({"b": pt_oe},
                                     deter_vals=obs.obsval.to_dict(),
                                     bins=40,
                                     filename=os.path.join(m_d, "pulu_glm_obs_summary.pdf"))
    pyemu.plot_utils.ensemble_res_1to1(pst=pst, ensemble={"b": pt_oe},
                                       filename=os.path.join(m_d, "pulu_glm_obs_vs_sim.pdf"))
    plt.show()


if __name__ == "__main__":
    start=time()
    setup_pulu()
    pulu_run_prior_mc()
    end=time()

    # start=time()
    # pulu_run_ies()
    # pulu_plot_ies_results()
    # end=time()
 
    start=time()
    pulu_run_glm()
    pulu_plot_glm_results()
    end=time()  
    print("total execution=",end-start)