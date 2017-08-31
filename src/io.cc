/* --------------------------------------------------------------------

   I/O routines that read the input file and HDF5 related routines
   Written by J. de la Cruz Rodriguez (ISP-SU 2016)

   Depends on HDF5 C routines with parallel MPI I/O installed

   -------------------------------------------------------------------- */

#include "hdf5.h"
#include "hdf5_hl.h"

#include <cstdio>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <mpi.h>
#include <dirent.h>
#include <sys/stat.h>

#include "model.h"
#include "io.h"
#include "physical_consts.h"
//
using namespace std;
using namespace modl;

/* ---------------------------------------------------------------------------- */

#define Q_WING              20.0

/* ---------------------------------------------------------------------------- */

int getAnum(char *elem)
{
  // Returns the atomic number assigned to an element label (e.g., "Fe").
  // anum(H) = 1, ... and so on.
  
  int anum = 0;
  for(int ii=0; ii<99; ii++){
    if(!strcmp(elem, phyc::ELEM[ii])){
      break;
    }
    anum = ii+1;
  }
  return anum;
}

/* ---------------------------------------------------------------------------- */

double air2vac(const double lambda_air){
  // Init values
  double lambda_vacuum;
  if(lambda_air > 2000) lambda_vacuum = lambda_air / 1.00029;
  else lambda_vacuum = lambda_air;

  // Iterate
  double error = 1.0;
  while(error > 1.e-6){
    error = lambda_air - vac2air(lambda_vacuum);
    lambda_vacuum = lambda_vacuum + error / 1.0029;
  }
  return lambda_vacuum;
}
/* ---------------------------------------------------------------------------- */

double vac2air(const double alamb){
  if(alamb < 2000.) return alamb;
  else return alamb/(1.0+2.735182e-4+131.4182/alamb/alamb+ 
		     2.76249e8/alamb/alamb/alamb/alamb);
}

/* ---------------------------------------------------------------------------- */

static const std::string vnames[13] = {
  "ltau500", "z", "temperature", "vz",
  "vx", "vy", "vturb", "B_str", "B_inc",
  "B_azi", "Pgas", "dens", "xne"};

/* ---------------------------------------------------------------------------- */

std::string removeSpaces(std::string input){
  input.erase(std::remove(input.begin(),input.end(),' '),input.end());
  return input;
}

/* ---------------------------------------------------------------------------- */

std::vector<std::string> strsplit(std::string &var, std::string token, bool rmspaces){
  
  std::vector<std::string> res;
  std::stringstream ss(var);

  std::string tmp;
  if(rmspaces)while(std::getline(ss, tmp, *(token.c_str()))) res.push_back(removeSpaces(tmp));
  else while(std::getline(ss, tmp, *(token.c_str()))) res.push_back(tmp);
  
  return res;
}

/* ---------------------------------------------------------------------------- */

void readInput(string filename, info &input, bool verbose)
{

  input.mu = 1.0, input.verbose = 1, input.nx=0, input.ny = 0, input.ndep = 0,
    input.nw=0, input.nstokes = 0, input.solver = 0, input.ipix = 0, input.log = stderr,
    input.temperature_cut = -1.0, input.eos_type = 0, input.gravity = 4.44,
    input.dlam = 2.0;
  
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in){
    std::string iline;
    while (std::getline(in, iline)) {
      // Split string at "=" and remove empty spaces
      std::string::size_type n;
      n = iline.find("=");
      std::string key = removeSpaces(iline.substr(0, n));
      std::string field = removeSpaces(iline.substr(n+1));
      bool set = false;
      
      // Fill in input structure
      if     (key == "input_model"){
	input.m.filename = field;
	set = true;
      }
      else if(key == "output_profiles"){
	input.p.filename = field;
	set = true;
      }
      else if(key == "logfile"){
	std::string lname = "log/logProc_"+std::to_string(input.myrank)+".log";

	if     (field == "stderr") input.log = stderr;
	else if(field == "stdout") input.log = stdout;
	else if(field == "master")
	  input.log = (input.myrank >0) ? fopen(lname.c_str(), "w") : stderr;
	else if(field == "null") input.log = fopen("/dev/null", "w");
	else{
	  if(!bdir_exists("log")) mkdir("log", 0700);
	  input.log = fopen(lname.c_str(), "w");
	  set = true;
	}
      }
      else if(key == "rt_solver"){
	input.solver = atoi(field.c_str());
	set = true;
      }
      else if(key == "eos"){
	if(field == std::string("pisk"))
	  input.eos_type = 1;
	set = true;
      }
      else if(key == "mu"){
	input.mu = atof(field.c_str());
	set = true;
      }
      else if(key == "gravity"){
	input.gravity = atof(field.c_str());
	set = true;
      }
      else if(key == "lines_file"){
	input.lines_file = field;
	set = true;
      }
      else if(key == "verbose"){
	input.verbose = atoi(field.c_str());
	set = true;
      }
      else if(key == "temperature_cut"){
	input.temperature_cut = atoi(field.c_str());
	set = true;
      }
      else if(key == "line_window"){
	input.dlam = atof(field.c_str());
	set = true;
      }
      else if(key == "region"){
	std::vector<std::string> param = strsplit(field,",");
	region reg = {};
	reg.w0 = air2vac(std::stod(param[0])); // convert to vacuum;
	reg.dw = std::stod(param[1]);
	reg.nw = std::stoi(param[2]);
	reg.cscal = std::stod(param[3]);
	//
	if(param.size() == 6){
	  reg.inst = param[4];
	  reg.ifile = param[5];
	}else{
	  reg.inst = "none";
	  reg.ifile = "none";
	}
  	input.nw += reg.nw;
	input.reg.push_back(reg);
	set = true;
      }
      else if(key == "") set = false;
      else{
	//	if(verbose) std::cout << "info: readInput: ignoring line: " << iline<<std::endl;
      }
      if(set && verbose) std::cout << "info: readInput: " << key << " -> "<<field << std::endl;
    }
  }

  in.close();
}

/* ---------------------------------------------------------------------------- */

bool bfile_exists(const std::string& name) {
  std::ifstream f(name.c_str(),std::ifstream::in);
  f.close();
  return f.good() ? true : false;
}

/* ---------------------------------------------------------------------------- */

bool bdir_exists(std::string name){
  DIR* dir = opendir(name.c_str());
  if (dir) return true;
  else if (ENOENT == errno) return false;
  else{
    fprintf(stderr,"error: bdir_exists: cannot write to folder [%s], exiting.\n", name.c_str());
    MPI_Abort(MPI_COMM_WORLD, 2);
    return true; // Just ensure function will return a value
  }
}

/* ---------------------------------------------------------------------------- */

bool h5varExists(hid_t &fid, std::string vname, hid_t &vid, FILE *log)
{
  bool res = (bool)H5LTfind_dataset(fid, vname.c_str());
  if(res){
    vid = H5Dopen2(fid, vname.c_str(), H5P_DEFAULT);
    fprintf(log, "[ %s ]",vname.c_str());
  }
  return res;
}

/* ---------------------------------------------------------------------------- */

void initReadIO(info &inp)
{

  const char routineName[] = "initReadIO";
  hid_t plist_id;

  
  /* --- Does input model exist ? --- */

  if(!bfile_exists(inp.m.filename)){
    fprintf(stderr,"error: initReadIO: file [%s] does not exist, exiting \n",
	    inp.m.filename.c_str());
    MPI_Abort(inp.comm, 2);
  }

  
  /* --- open file, setting parallel I/O flags --- */
  
  if ((plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0)
    h5err(inp.myrank, routineName);

  if ((H5Pset_fapl_mpio(plist_id, inp.comm, inp.info)) < 0)
    h5err(inp.myrank, routineName);

  if ((inp.m.fid = H5Fopen(inp.m.filename.c_str(), H5F_ACC_RDONLY, plist_id)) < 0)
    h5err(inp.myrank, routineName);

  if ((H5Pclose(plist_id)) < 0)
    h5err(inp.myrank, routineName); 


  
  /* --- Get dims --- */
  
  int ndim = 0;
  if((H5LTget_dataset_ndims(inp.m.fid, "temperature",  &ndim) < 0))
    h5err(inp.myrank, routineName);
  
  hsize_t  dims[ndim];
  if ((H5LTget_dataset_info(inp.m.fid, "temperature", dims, NULL, NULL)) < 0)
    h5err(inp.myrank, routineName);


  if(ndim != 4){
    fprintf(inp.log, "error: initReadIO: the model variables must have 4 dimensions [nt, ny, nx, ndep], exiting.\n");
    h5err(inp.myrank, routineName); 
  }else{
    inp.nt = dims[0], inp.ny = dims[1], inp.nx = dims[2], inp.ndep = dims[3];
  }
  
  
  if(inp.verbose){
    fprintf(inp.log,"info: initReadIO: Found dims [ nt=%d, ny=%d, nx=%d, ndep=%d ]\n", inp.nt, inp.ny, inp.nx, inp.ndep);
  }

  
  /* --- Init variable IDs --- */
  
  if((H5LTget_attribute_int(inp.m.fid, "/", "units", &inp.units)) < 0) inp.units = 0;
  inp.vdef.resize(13,false);

  fprintf(inp.log, "info: initReadIO: Found vars ");
  for(size_t ii=0; ii<13; ii++)
    inp.vdef[ii] =  h5varExists(inp.m.fid, vnames[ii], inp.m.vid[ii],  inp.log);
  fprintf(inp.log, "\n");

  
  /* --- Init mem space and dataspace --- */

  hsize_t dmem = (hsize_t)inp.ndep;
  if ((inp.m.mid = H5Screate_simple(1, &dmem, NULL)) < 0) h5err(inp.myrank, routineName);
  if ((inp.m.did = H5Dget_space(inp.m.vid[2])) < 0) h5err(inp.myrank, routineName);
  
}

/* ---------------------------------------------------------------------------- */

void h5err(int myrank, const char *rname)
{  
  fprintf(stderr,"error: Process %3d: %s: HDF5 error.\n", myrank, rname);
  MPI_Abort(MPI_COMM_WORLD, 2);
}

/* ---------------------------------------------------------------------------- */

void readPix(hid_t &hid, hsize_t msp, hsize_t  dsp, double *var, FILE *log)
{
  if ((H5Dread(hid, H5T_NATIVE_DOUBLE, msp, dsp, H5P_DEFAULT, var)) < 0){
    fprintf(log, "error: readPix: cannot read variable, exiting!\n");
    MPI_Abort(MPI_COMM_WORLD, 2);
  }else return;
}

  
/* ---------------------------------------------------------------------------- */

void readAtmosTYX(size_t tt, size_t yy, size_t xx, mdepth &m, info &inp)
{

  const char routineName[] = "readAtmosTYX";
  
  /* --- Check dims and re-allocate if needed --- */

  if(m.ndep != inp.ndep){
    m.init(inp.ndep, true, (modl::munit)inp.units);
  }
  

  /* --- compute file offset and read vars --- */

  hsize_t start[]    = {tt, yy, xx, 0}; 
  hsize_t count[]    = {1, 1, 1, (hsize_t)inp.ndep};


  /* --- Select pixel to read --- */

  if ((H5Sselect_hyperslab(inp.m.did, H5S_SELECT_SET, start, NULL, count, NULL)) < 0) h5err(inp.myrank, routineName); 


  
  /* --- Read defined variables in the model --- */
  
  for(size_t ii=0; ii<13; ii++){
    if(inp.vdef[ii])
      readPix(inp.m.vid[ii], inp.m.mid, inp.m.did, &m.buf(ii,0), inp.log);
  }
  if(inp.verbose >= 2) fprintf(inp.log, "info: read nt=%4zu, ny=%4zu, nx=%4zu\n", tt, yy, xx);
  
}

/* ---------------------------------------------------------------------------- */

void writeProfileTYX(size_t tt, size_t yy, size_t xx, double *sp, info &inp)
{
  hsize_t    offset[4] = {0,0,0,0}, count[4] = {1,1,1,1};
  count[3] = (hsize_t)inp.nw;
  offset[0] = (hsize_t)tt, offset[1] = (hsize_t)yy, offset[2] = (hsize_t)xx;
  
  const char routineName[] = "writeProfileTYX";

  /* --- Select hyperslab in file --- */
  
  if (( H5Sselect_hyperslab(inp.p.did, H5S_SELECT_SET, offset, NULL, count, NULL) ) < 0)
    h5err(inp.myrank, routineName);


  /* --- Write data --- */
  
  if (( H5Dwrite(inp.p.vid[0], H5T_NATIVE_DOUBLE, inp.p.mid, inp.p.did, H5P_DEFAULT, sp) ) < 0) 
    h5err(inp.myrank, routineName); 
}

/* ---------------------------------------------------------------------------- */

void initWriteIO(info &inp, double *lambda)
{
  hid_t plist;
  const char routineName[] = "initWriteIO";


  /* --- Open file and set attributes --- */
  
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0)  h5err(inp.myrank, routineName); 
  if (( H5Pset_fapl_mpio(plist, inp.comm, inp.info) ) < 0)  h5err(inp.myrank, routineName); 
  if (( inp.p.fid = H5Fcreate(inp.p.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist) ) < 0)
    h5err(inp.myrank, routineName); 
  if (( H5Pclose(plist) ) < 0)  h5err(inp.myrank, routineName); 

  

  /* --- Create dataspace and datasets, data are stored as float32 --- */
  
  hsize_t dims[4];
  dims[0] = (hsize_t)inp.nt, dims[1] = (hsize_t)inp.ny, dims[2] = (hsize_t)inp.nx,
    dims[3] = (hsize_t)inp.nw;
  
  if (( inp.p.did = H5Screate_simple(4, dims, NULL) ) < 0) h5err(inp.myrank, routineName);
  if (( plist = H5Pcreate(H5P_DATASET_CREATE) ) < 0) h5err(inp.myrank, routineName);
  if (( inp.p.vid[0] = H5Dcreate(inp.p.fid, "Stokes_I", H5T_NATIVE_FLOAT,inp.p.did,
				 H5P_DEFAULT, plist, H5P_DEFAULT)) < 0)
    h5err(inp.myrank, routineName);


  
  /* --- Store wavelength array --- */
  
  dims[0] = inp.nw;
  if (( H5LTmake_dataset(inp.p.fid, "Wavelength", 1, dims, H5T_NATIVE_DOUBLE, lambda) ) < 0)
    h5err(inp.myrank, routineName);
  
  if (( H5LTset_attribute_string(inp.p.fid, "Wavelength", "units",
                                 "Angstroms") ) < 0) h5err(inp.myrank, routineName);
  


  /* --- Init memspace and dataspace --- */

  dims[0] = (hsize_t)inp.nw;
  if (( inp.p.mid = H5Screate_simple(1, dims, NULL) ) < 0) h5err(inp.myrank, routineName);
  if (( inp.p.did = H5Dget_space(inp.p.vid[0]) ) < 0)  h5err(inp.myrank, routineName);
  
  if (( H5Pclose(plist) ) < 0)  h5err(inp.myrank, routineName); 

}

/* ---------------------------------------------------------------------------- */

void closeProfile(info &inp)
{

  const char routineName[] = "closeProfile";
  
  if (( H5Sclose(inp.p.did) ) < 0) h5err(inp.myrank, routineName); 
  if (( H5Sclose(inp.p.mid) ) < 0) h5err(inp.myrank, routineName);
  //
  if (( H5Dclose(inp.p.vid[0]) ) < 0) h5err(inp.myrank, routineName); 
  if (( H5Fclose(inp.p.fid) ) < 0) h5err(inp.myrank, routineName); 
}

/* ---------------------------------------------------------------------------- */

void closeModel(info &inp)
{
  const char routineName[] = "closeModel";
    
  if (( H5Sclose(inp.m.did) ) < 0)  h5err(inp.myrank, routineName); 
  if (( H5Sclose(inp.m.mid) ) < 0)  h5err(inp.myrank, routineName); 

  for(size_t ii=0; ii<13; ii++) if(inp.vdef[ii]) H5Dclose(inp.m.vid[ii]);
  H5Fclose(inp.m.fid);
}

/* ---------------------------------------------------------------------------- */

int readValdLines(string filename, info &input)
{
  
  /* --- Open file --- */

  int ncounted = 0, nlines = 0;
  bool firsttime = true;
  
  std::ifstream in(filename, std::ios::in | std::ios::binary);

  if (in){
    std::string iline;
    while (std::getline(in, iline)) {
      if(iline.size() < 1) continue;
      if(iline[0] == '#' || iline[0] == '=' || iline[0] == 'W') continue;

      if(firsttime){ // Get number of lines
	std::vector<std::string> dum  = strsplit(iline, ",", true);
	nlines = std::atoi(dum[2].c_str());
	firsttime = false;
	break;
      }
    }

    /* --- Now we know the number of lines, and there are 4 lines in the file 
       per spectral line --- */

    for(int ii = 0; ii<nlines; ii++){
      line_t iline = {};
      
      for(int ll = 0; ll<4; ll++){
	std::string tmp;
	while(std::getline(in, tmp) && tmp[0] != '\''); // Get next valid line

	std::vector<std::string> els = strsplit(tmp, ",", true);
	
	switch(ll){
	case 0:
	  {
	    // Element name
	    strcpy(iline.elem, els[0].substr(1,2).c_str());
	    iline.anum = getAnum(iline.elem);
	    iline.amass = phyc::AMASS[iline.anum];
	    // Ion stage
	    iline.ion = (int)std::atoi(els[0].substr(3,1).c_str());
	    // lambda_0 (air)
	    iline.w0 = air2vac((double)std::atof(els[1].c_str()));
	    iline.nu0 = (phyc::CC * 1.0E8)  / iline.w0;  // in Hz (s^-1)
	    iline.width = input.dlam;
	    // gf
	    iline.gf = pow(10., (double)std::atof(els[2].c_str()));
	    // E_low
	    iline.e_low = (double)std::atof(els[3].c_str());
	    // j_low
	    iline.Jlow = (double)std::atof(els[4].c_str());
	    // E_up
	    iline.e_up = (double)std::atof(els[5].c_str());
	    // j_up
	    iline.Jup = (double)std::atof(els[6].c_str());
	    // G_low (Lande)
	    iline.Glow = (double)std::atof(els[7].c_str());
	    // G_up (Lande)
	    iline.Gup = (double)std::atof(els[8].c_str());
	    // Do not read mean, not needed
	    // Rad. damp const.
	    iline.g_rad = (double)std::atof(els[10].c_str());
	    // Stark damp. const.
	    iline.g_str = (double)std::atof(els[11].c_str());
	    // VdW damp. const.
	    iline.g_vdw = (double)std::atof(els[12].c_str());
	    
	    /* --- Manipulate units and fill other quantities --- */
	    if(iline.ion == 1) iline.eion = phyc::EION1[iline.anum] - iline.e_low;
	    else if(iline.ion == 2) iline.eion = phyc::EION2[iline.anum] - iline.e_low;
	    else iline.eion = 0.0;

	    iline.e_up  *= phyc::EV;
	    iline.e_low *= phyc::EV;
	    iline.eion  *= phyc::EV;

	    // Fix gamma rad
	    if(iline.g_rad == 0.0) iline.g_rad = 0.22 / (iline.w0*iline.w0*1.e-16); // Gray (2005) eq. 11.13
	    else iline.g_rad = pow(10.0,iline.g_rad);

	    // gamma Stark (Vald gives log10(gamma_4/nne at 10000. [K]).
	    if(iline.g_str != 0.0) iline.g_str = pow(10.0, iline.g_str);
	    else iline.g_str = 0.0;

	    // vdW (Barklem or Unsold)
	    if(iline.g_vdw > 20.0){ // Barklem
	      iline.b_sig = int(iline.g_vdw);
	      iline.b_alp = iline.g_vdw - int(iline.g_vdw);
	      iline.barklem = true;
	    } else{
	      iline.b_sig = 0.0;
	      iline.b_alp = 0.0;
	      iline.barklem = false;
	      iline.g_vdw = pow(10.0, iline.g_vdw);
	    }
	    iline.firsttime = true;
	    input.lin.push_back(iline);

	    break;
	  }
	default:
	  break;
	} //case
      } // ll
    } // ss
  } // in

  fprintf(input.log, "info: io, read_lines: read [%d] lines\n", (int)input.lin.size());
}
