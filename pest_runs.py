import os
from time import time
import shutil
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
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
    output_filename = os.path.join(case, "output_climax_misti.txt")
    ins_filename, obs_df = write_output_ins_file(output_filename)
    config_filename = os.path.join(os.path.join(case, "config.py"))
    tpl_filename, par_df = write_config_tpl_file(config_filename)

    # create a control file instance
    pst = pyemu.Pst.from_io_files(tpl_files=tpl_filename, in_files=config_filename, ins_files=ins_filename,
                                  out_files=output_filename, pst_path=".")
    pst.model_input_data.iloc[0,:].loc["model_file"] = os.path.join("src",os.path.split(config_filename)[-1])
    # set the observation values
    pst.observation_data.loc[obs_df.obsnme, "obsval"] = obs_df.obsval

    # set the parameter information for the pulu pars
    par = pst.parameter_data
    # start by fixing all the parameters in the control file
    par.loc[:, "partrans"] = "fixed"

    par.loc["total_erupted_mass", "parlbnd"] = 0.1e11
    par.loc["total_erupted_mass", "parubnd"] = 5.0e11
    par.loc["total_erupted_mass", "parval1"] = 1.0e11
    par.loc["total_erupted_mass", "partrans"] = "none"

    par.loc["column_height", "parlbnd"] = 10000
    par.loc["column_height", "parubnd"] = 30000
    par.loc["column_height", "parval1"] = 25000
    par.loc["column_height", "partrans"] = "none"

    par.loc["diffusion_coef", "parlbnd"] = 300
    par.loc["diffusion_coef", "parubnd"] = 3000
    par.loc["diffusion_coef", "parval1"] = 1500
    par.loc["diffusion_coef", "partrans"] = "none"

    # it looks like mean grainsize is already in log space?
    par.loc["tgsd_mean", "parlbnd"] = -2
    par.loc["tgsd_mean", "parubnd"] = 2
    par.loc["tgsd_mean", "parval1"] = -1.99
    par.loc["tgsd_mean", "partrans"] = "none"
    par.loc["tgsd_mean", "parchglim"] = "relative"

    par.loc["ellipse_major_axis", "parlbnd"] = 5000
    par.loc["ellipse_major_axis", "parubnd"] = 25000
    par.loc["ellipse_major_axis", "parval1"] = 11000
    par.loc["ellipse_minor_axis", "parval1"] = 11000

 
    # tie the minor axis par to the major so they function
    # as a single radius
    par.loc["ellipse_major_axis", "partrans"] = "none"
    par.loc["ellipse_minor_axis", "partrans"] = "tied"
    par.loc["ellipse_minor_axis", "partied"] = "ellipse_major_axis"
  
    par.loc["wind_speed", "parlbnd"] = 0.5
    par.loc["wind_speed", "parubnd"] = 20.0
    par.loc["wind_speed", "parval1"] = 5.0
    par.loc["wind_speed", "partrans"] = "none"

    par.loc["wind_direction", "parlbnd"] = 160
    par.loc["wind_direction", "parubnd"] = 260
    par.loc["wind_direction", "parval1"] = 220
    par.loc["wind_direction", "partrans"] = "none"

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


def run_glm(case,noptmax=5,num_reals=3000,num_workers=10):
    """run pestpp-glm in parallel locally"""
    pst = pyemu.Pst(os.path.join(case, case+".pst"))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["glm_num_reals"] = num_reals

    pst.write(os.path.join(case, case+"_run.pst"))
    pyemu.os_utils.start_workers(case, "pestpp-glm", case+"_run.pst", num_workers=num_workers,
                                 master_dir=case+"_glm_master")


def run_prior_monte_carlo(case,num_reals=3000,num_workers=10,use_uniform=True):
    pst = pyemu.Pst(os.path.join(case, case+".pst"))

    if use_uniform:
        pst.parameter_data.loc[pst.adj_par_names,"partrans"] = "none"
        pe = pyemu.ParameterEnsemble.from_uniform_draw(pst,num_reals=num_reals)
        pe.to_csv(os.path.join(case,case+"_prior.csv"))
        pst.pestpp_options["ies_par_en"] = case+"_prior.csv"

    pst.control_data.noptmax = -1
    pst.pestpp_options["ies_num_reals"] = num_reals

    pst.write(os.path.join(case, case+"_run.pst"))
    pyemu.os_utils.start_workers(case, "pestpp-ies", case+"_run.pst", num_workers=num_workers,
                                 master_dir=case+"_pmc_master")

def plot_glm_results(case,pmc_dir=None):
    """plot the pestpp-glm results"""
    m_d = case+"_glm_master"
    #print('Done')
    pst = pyemu.Pst(os.path.join(m_d, case+"_run.pst"))
    #print(pst.phi)
    
    if pmc_dir is not None:
        pr_oe = pd.read_csv(os.path.join(pmc_dir, case+"_run.0.obs.csv"), index_col=0)
        pr_oe = pyemu.ObservationEnsemble(pst=pst,df=pr_oe)
        pr_pe = pd.read_csv(os.path.join(pmc_dir, case+"_run.0.par.csv"), index_col=0)
        pr_pe = pr_pe.loc[:,pst.adj_par_names]


    pt_pe = pd.read_csv(os.path.join(m_d, case+"_run.post.paren.csv"), index_col=0)
    pt_oe = pd.read_csv(os.path.join(m_d, case+"_run.post.obsen.csv"), index_col=0)
    pt_oe = pyemu.ObservationEnsemble(pst=pst,df=pt_oe)

    # rejection sampling - only keep posterior realizations that are within XXX% of the best phi
    pt_pv = pt_oe.phi_vector
    best_phi = min(pst.phi,pt_pv.min())
    acc_phi = best_phi * 1.1
    pt_pv = pt_pv.loc[pt_pv<acc_phi]
    pt_pe = pt_pe.loc[pt_pv.index,:]
    pt_oe = pt_oe.loc[pt_pv.index,:]
    print("glm phi:",pst.phi, "best phi:",best_phi,"passing realizations:",pt_oe.shape[0])
    fig,ax = plt.subplots(1,1,figsize=(4,4))
    if pmc_dir is not None:
        ax.hist(pr_oe.phi_vector.apply(np.log10),bins=20,facecolor="0.5",alpha=0.5,edgecolor="none",density=False)
    ax.hist(pt_oe.phi_vector.apply(np.log10),bins=20,facecolor="b",alpha=0.5,edgecolor="none",density=False)
    ax.plot([pst.phi,pst.phi],ax.get_ylim(),"b--")
    ax.set_title("best phi:{0:5.2E}, acceptable phi:{1:5.2E}, number of realizations passing: {2}".format(best_phi,acc_phi,pt_oe.shape[0]))
    ax.set_xlabel("$log_{10} \phi$")
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
                                         filename=os.path.join(m_d, case+"_glm_obs_summary.pdf"),
                                         std_window=0.5,deter_range=True,sync_bins=False)
        pyemu.plot_utils.ensemble_res_1to1(pst=pst, ensemble={"0.5":pr_oe,"b": pt_oe},alpha=0.5,
                                       filename=os.path.join(m_d, case+"_glm_obs_vs_sim.pdf"),
                                       base_ensemble=os.path.join(pmc_dir,case+"_run.obs+noise.csv"))
        pst.parameter_data.loc[:,"pargp"] = pst.par_names
        pyemu.plot_utils.ensemble_change_summary(pst=pst, ensemble1=pr_pe,ensemble2=pt_pe,
                                       filename=os.path.join(m_d, case+"_glm_par_change_summary.pdf"))
    #plt.show()


def sensitivity_experiment():
    case = "misti"
    pst = pyemu.Pst(os.path.join(case, case + ".pst"))
    pst.control_data.noptmax = 0
    pst.write(os.path.join(case, case + ".pst"))
    pyemu.os_utils.run("pestpp-glm {0}.pst".format(case),cwd=case)
    pst = pyemu.Pst(os.path.join(case, case + ".pst"))
    base_phi = pst.phi
    pst.parameter_data.loc["wind_speed","parval1"] *= 100
    pst.write(os.path.join(case, case + ".pst"))
    pyemu.os_utils.run("pestpp-glm {0}.pst".format(case), cwd=case)
    pst = pyemu.Pst(os.path.join(case, case + ".pst"))
    pert_phi = pst.phi
    print("base phi:",base_phi,"perturbed phi:", pert_phi)


def plot_glue_results(case,glm_dir=None):
    """plot the rejection sampling results.  If glm_dir is not None, then glm-based fosm
    distributions are also plotted

    """

    m_d = case+"_pmc_master"
    pst = pyemu.Pst(os.path.join(m_d, case+"_run.pst"))

    fosm_df = None
    glm_pst = None
    if glm_dir is not None:
        glm_pst = pyemu.Pst(os.path.join(glm_dir, case + "_run.pst"))
        par = glm_pst.parameter_data
        log_pars = par.loc[par.partrans == "log", "parnme"].values
        assert os.path.exists(glm_dir)
        fosm_punc_file = os.path.join(glm_dir,case+"_run.par.usum.csv")
        assert os.path.exists(fosm_punc_file),fosm_punc_file
        fosm_df = pd.read_csv(fosm_punc_file,index_col=0)
        fosm_df.index = fosm_df.index.map(str.lower)
        fosm_df.loc[log_pars,:] = 10.0**fosm_df.loc[log_pars,:]

       
    pr_oe = pd.read_csv(os.path.join(m_d, case+"_run.0.obs.csv"), index_col=0)

    pr_pe = pd.read_csv(os.path.join(m_d, case+"_run.0.par.csv"), index_col=0)
    pr_pe.index = pr_pe.index.map(str)
    pr_oe.index = pr_oe.index.map(str)
    pr_pe = pr_pe.loc[pr_oe.index,:]
    pr_oe = pyemu.ObservationEnsemble(pst=pst, df=pr_oe)
    pr_pe = pr_pe.loc[:,pst.adj_par_names]

    print(pr_oe.shape,pr_pe.shape)

  # rejection sampling - only keep posterior realizations that are within XXX% of the best phi
    pr_pv = pr_oe.phi_vector
    best_phi = min(pst.phi,pr_pv.min())
    acc_phi = best_phi * 2
    pt_pv = pr_pv.loc[pr_pv<acc_phi]
    pt_oe = pr_oe.loc[pt_pv.index,:]
    pt_pe = pr_pe.loc[pt_oe.index,:]
    print("best phi:",best_phi,"passing realizations:",pt_oe.shape[0])
    pt_pe.to_csv("filtered_posterior.csv")
    pr_pe.to_csv("prior.csv")

    fig,ax = plt.subplots(1,1,figsize=(4,4))
    
    ax.hist(pr_oe.phi_vector.apply(np.log10),bins=100,facecolor="0.5",alpha=0.5,edgecolor="none",density=False,label="prior")
    ax.hist(pt_oe.phi_vector.apply(np.log10),bins=100,facecolor="b",alpha=0.5,edgecolor="none",density=False,label="posterior")
    ylim = ax.get_ylim()
    if glm_pst is not None:
        ax.plot([np.log10(glm_pst.phi), np.log10(glm_pst.phi)], ylim, "r--",label="glm phi")
    ax.plot([np.log10(best_phi),np.log10(best_phi)],ylim,"b--",label="GLUE best phi")
    ax.set_ylim(ylim)
    ax.set_title("best phi:{0:5.2E}, acceptable phi:{1:5.2E}, number of realizations passing: {2}".format(best_phi,acc_phi,pt_oe.shape[0]))
    ax.set_xlabel("$log_{10} \phi$")
    ax.legend(loc="upper right")
    plt.savefig(os.path.join(m_d,"phi_hist.pdf"))
    plt.close(fig)
  

    #pyemu.plot_utils.ensemble_helper({"b": pt_pe,"0.5":pr_pe}, bins=100,
    #                                 filename=os.path.join(m_d, case+"_summary.pdf"))
    #plt.show()
    bins = 20
    with PdfPages(os.path.join(m_d,"par_summary.pdf")) as pdf:
        for par_name in pst.adj_par_names:
            fig,ax = plt.subplots(1,1,figsize=(6,3))
            ax.hist(pr_pe.loc[:,par_name],bins=bins,alpha=0.5,edgecolor="none",facecolor="0.5",density=True,label="prior")
            ax.hist(pt_pe.loc[:, par_name], bins=bins, alpha=0.5, edgecolor="none", facecolor="b",density=True,label="posterior")
            ax.set_title(par_name,loc="left")
            if fosm_df is not None:
                axt = plt.twinx(ax)
                print(par_name,fosm_df.loc[par_name,"prior_mean"],fosm_df.loc[par_name,"prior_stdev"])
                pr_x,pr_y = pyemu.plot_utils.gaussian_distribution(fosm_df.loc[par_name,"prior_mean"],fosm_df.loc[par_name,"prior_stdev"])
                axt.plot(pr_x,pr_y,color="0.5",dashes=(1,1),lw=3,label="FOSM prior")
                pt_x, pt_y = pyemu.plot_utils.gaussian_distribution(fosm_df.loc[par_name, "post_mean"],
                                                                    fosm_df.loc[par_name, "post_stdev"])
                axt.plot(pt_x, pt_y, color="b", dashes=(1, 1), lw=3,label="FOSM posterior")
                axt.set_yticks([])
                ylim = ax.get_ylim()
                xlim = ax.get_xlim()
                ax.plot([0,0],[0,0],color="b", dashes=(1, 1), lw=3,label="FOSM posterior")
                ax.plot([0, 0], [0, 0], color="0.5", dashes=(1, 1), lw=3, label="FOSM prior")
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)

            ax.legend(loc="upper right")
            plt.tight_layout()
            pdf.savefig()
            plt.close(fig)

    obs = pst.observation_data
    
    pyemu.plot_utils.ensemble_helper({"b": pt_oe,"0.5":pr_oe},
                                     deter_vals=obs.obsval.to_dict(),
                                     bins=100,
                                     filename=os.path.join(m_d, case+"_glm_obs_summary.pdf"),
                                     std_window=0.5,deter_range=True,sync_bins=False)
    pyemu.plot_utils.ensemble_res_1to1(pst=pst, ensemble={"0.5":pr_oe,"b": pt_oe},alpha=0.5,
                                   filename=os.path.join(m_d, case+"_glm_obs_vs_sim.pdf"),
                                   base_ensemble=os.path.join(m_d,case+"_run.obs+noise.csv"))
    pst.parameter_data.loc[:,"pargp"] = pst.par_names
    pyemu.plot_utils.ensemble_change_summary(pst=pst, ensemble1=pr_pe,ensemble2=pt_pe,
                                   filename=os.path.join(m_d, case+"_glm_par_change_summary.pdf"))
    #plt.show()


if __name__ == "__main__":
    volcano = "misti" # working directory with volcano data

    start=time()
    #setup(volcano)
    #sensitivity_experiment()
    #run_prior_monte_carlo(volcano,num_reals=500,num_workers=10)
    run_glm(volcano,noptmax=10,num_reals=0,num_workers=15)
    #plot_glm_results(volcano,pmc_dir="{0}_pmc_master".format(volcano))
    plot_glue_results(volcano,glm_dir="misti_glm_master")
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