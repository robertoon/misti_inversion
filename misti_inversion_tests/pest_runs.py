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


def setup(case):
    # run the model one
    shutil.copy2(os.path.join(case, "org_config.py"), os.path.join(case, "config.py"))
    pyemu.os_utils.run("python main.py", cwd=case)

    # prepare the interface files
    output_filename = os.path.join(case, "output_misti.txt")
    ins_filename, obs_df = write_output_ins_file(output_filename)
    config_filename = os.path.join(os.path.join(case, "config.py"))
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

    par.loc["total_erupted_mass", "parlbnd"] = 1.3e10
    par.loc["total_erupted_mass", "parubnd"] = 6.8e10
    par.loc["total_erupted_mass", "parval1"] = 6.3e10
    par.loc["total_erupted_mass", "partrans"] = "log"

    par.loc["column_height", "parlbnd"] = 21000
    par.loc["column_height", "parubnd"] = 31000
    par.loc["column_height", "parval1"] = 26000
    par.loc["column_height", "partrans"] = "log"

    par.loc["diffusion_coef", "parlbnd"] = 1000
    par.loc["diffusion_coef", "parubnd"] = 2000
    par.loc["diffusion_coef", "parval1"] = 1500
    par.loc["diffusion_coef", "partrans"] = "log"

    # it looks like mean grainsize is already in log space?
    par.loc["tgsd_mean", "parlbnd"] = -6
    par.loc["tgsd_mean", "parubnd"] = 6
    par.loc["tgsd_mean", "parval1"] = -1.99
    par.loc["tgsd_mean", "partrans"] = "none"
    par.loc["tgsd_mean", "parchglim"] = "relative"

    par.loc["ellipse_major_axis", "parlbnd"] = 5000
    par.loc["ellipse_major_axis", "parubnd"] = 18000
    par.loc["ellipse_major_axis", "parval1"] = 11500
    par.loc["ellipse_minor_axis", "parval1"] = 11500

    # tie the minor axis par to the major so they function
    # as a single radius
    par.loc["ellipse_major_axis", "partrans"] = "none"
    par.loc["ellipse_minor_axis", "partrans"] = "tied"
    par.loc["ellipse_minor_axis", "partied"] = "ellipse_major_axis"
  
    par.loc["wind_speed", "parlbnd"] = 1.0
    par.loc["wind_speed", "parubnd"] = 7.0
    par.loc["wind_speed", "parval1"] = 4.0
    par.loc["wind_speed", "partrans"] = "log"

    par.loc["wind_direction", "parlbnd"] = 210
    par.loc["wind_direction", "parubnd"] = 226
    par.loc["wind_direction", "parval1"] = 218
    par.loc["wind_direction", "partrans"] = "log"

    # set the observed thickness values
    pst.try_parse_name_metadata()  # this sets individual columns for northing and easting

    pst.write(os.path.join(case, case+".pst"))

    obs = pst.observation_data
    # read the observed data
    df = pd.read_csv(os.path.join(case, case+"_grainsize.csv"))
    df.columns = df.columns.map(lambda x: x.lower().replace(" ", "_"))

    # for each observation...
    for oname, oe, on in zip(obs.obsnme, obs.e, obs.n):
        # find the location based on easting and northing (not the relative, the absolute)
        print(oe,on)
        v = df.loc[df.apply(lambda x: x.easting == float(oe) and x.northing == float(on), axis=1), "thickness"].values[
            0]

        # set the observed value
        obs.loc[oname, "obsval"] = v

    # set the model forward run command
    pst.model_command = "python main.py"

    # run pestpp-glm for a single one-off forward run
    pst.control_data.noptmax = 0
    pst.write(os.path.join(case, case+".pst"))
    pyemu.os_utils.run("pestpp-glm {0}.pst".format(case), cwd=case)
    pst = pyemu.Pst(os.path.join(case, case+".pst"))
    pst.plot(kind="1to1")
    #plt.show()


def run_glm(case,noptmax=5,num_reals=3000):
    """run pestpp-glm in parallel locally"""
    pst = pyemu.Pst(os.path.join(case, case+".pst"))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["glm_num_reals"] = num_reals

    pst.write(os.path.join(case, case+"_run.pst"))
    pyemu.os_utils.start_workers(case, "pestpp-glm", case+"_run.pst", num_workers=10,
                                 master_dir=case+"_glm_master")


def run_prior_monte_carlo(case,num_reals=3000):
    pst = pyemu.Pst(os.path.join(case, case+".pst"))
    pst.control_data.noptmax = -1
    pst.pestpp_options["ies_num_reals"] = num_reals

    pst.write(os.path.join(case, case+"_run.pst"))
    pyemu.os_utils.start_workers(case, "pestpp-ies", case+"_run.pst", num_workers=10,
                                 master_dir=case+"_pmc_master")

def plot_glm_results(case,pmc_dir=None):
    """plot the pestpp-glm results"""
    m_d = case+"_glm_master"
    pst = pyemu.Pst(os.path.join(m_d, case+"_run.pst"))
    print(pst.phi)
    
    if pmc_dir is not None:
        pr_oe = pd.read_csv(os.path.join(pmc_dir, case+"_run.0.obs.csv"), index_col=0)
        pr_oe = pyemu.ObservationEnsemble(pst=pst,df=pr_oe)
        pr_pe = pd.read_csv(os.path.join(pmc_dir, case+"_run.0.par.csv"), index_col=0)
        pr_pe = pr_pe.loc[:,pst.adj_par_names]


    pt_pe = pd.read_csv(os.path.join(m_d, case+"_run.post.paren.csv"), index_col=0)
    pt_oe = pd.read_csv(os.path.join(m_d, case+"_run.post.obsen.csv"), index_col=0)
    pt_oe = pyemu.ObservationEnsemble(pst=pst,df=pt_oe)


    fig,ax = plt.subplots(1,1,figsize=(4,4))
    if pmc_dir is not None:
        ax.hist(pr_oe.phi_vector,bins=20,facecolor="0.5",alpha=0.5,edgecolor="none",density=False)
    ax.hist(pt_oe.phi_vector,bins=20,facecolor="b",alpha=0.5,edgecolor="none",density=False)
    ax.plot([pst.phi,pst.phi],ax.get_ylim(),"b--")
    plt.savefig(os.path.join(m_d,"phi_hist.pdf"))
    plt.close(fig)
  
    pt_pe = pt_pe.loc[:, pst.adj_par_names]

    if pmc_dir is None:
        pyemu.plot_utils.ensemble_helper({"b": pt_pe}, bins=100,
                                         filename=os.path.join(m_d, case+"_summary.pdf"))
    else:
        pyemu.plot_utils.ensemble_helper({"b": pt_pe,"0.5":pr_pe}, bins=100,
                                         filename=os.path.join(m_d, case+"_summary.pdf"))
    #plt.show()
    obs = pst.observation_data
    if pmc_dir is None:
        pyemu.plot_utils.ensemble_helper({"b": pt_oe},
                                         deter_vals=obs.obsval.to_dict(),
                                         bins=100,
                                         filename=os.path.join(m_d, case+"_glm_obs_summary.pdf"))
        pyemu.plot_utils.ensemble_res_1to1(pst=pst, ensemble={"b": pt_oe},
                                       filename=os.path.join(m_d, case+"_glm_obs_vs_sim.pdf"))
    else:
        pyemu.plot_utils.ensemble_helper({"b": pt_oe,"0.5":pr_oe},
                                         deter_vals=obs.obsval.to_dict(),
                                         bins=100,
                                         filename=os.path.join(m_d, case+"_glm_obs_summary.pdf"))
        pyemu.plot_utils.ensemble_res_1to1(pst=pst, ensemble={"0.5":pr_oe,"b": pt_oe},alpha=0.5,
                                       filename=os.path.join(m_d, case+"_glm_obs_vs_sim.pdf"),
                                       base_ensemble=os.path.join(pmc_dir,case+"_run.obs+noise.csv"))
    #plt.show()



if __name__ == "__main__":
    volcano = "misti" # working directory with volcano data

    start=time()
    #setup(volcano)
    run_prior_monte_carlo(volcano,num_reals=1000)
    run_glm(volcano,num_reals=1000)
    plot_glm_results(volcano,pmc_dir="{0}_pmc_master".format(volcano))
    end=time()
    print("total execution=",end-start)


    #start=time()
    #setup_pulu()
    #pulu_run_prior_mc()
    #end=time()

    # start=time()
    # pulu_run_ies()
    # pulu_plot_ies_results()
    # end=time()
 
    # start=time()
    #pulu_run_glm()
    
    # pulu_plot_glm_results()
    # end=time()  
    # print("total execution=",end-start)