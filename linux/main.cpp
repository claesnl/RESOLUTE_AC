/*
* Build RESOLUTE umap
*
* VERSION 	:: 	2.0.1
* BUILD 	::	26-July-2016
* AUTHOR	:: 	Claes Ladefoged
*
* Changes	::
*	v2.0.1	::	26-July-2016 :: Added automatic build scripts. Working build now for Ubuntu and CentOS.
*
*/

#include <iostream>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <sys/types.h>
#include <fcntl.h>
#include <dirent.h>
#include <vector>
#include <math.h>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "dcmtk/dcmdata/dcpxitem.h"

using namespace std;

#define WIDTH 192
#define HEIGHT 192
#define DEPTH 192

float *points;

float *ute1 = new float[ WIDTH*HEIGHT*DEPTH ];
float *ute2 = new float[ WIDTH*HEIGHT*DEPTH ];

float *sinusmask = new float[ WIDTH*HEIGHT*DEPTH ];
float *brain = new float[ WIDTH*HEIGHT*DEPTH ];
float *csf = new float[ WIDTH*HEIGHT*DEPTH ];
float *base = new float[ WIDTH*HEIGHT*DEPTH ];
float *r2noise = new float[ WIDTH*HEIGHT*DEPTH ];
float *sphenoidmask = new float[ WIDTH*HEIGHT*DEPTH ];

float *ute1_s = new float[ WIDTH*HEIGHT*DEPTH ];
float *ute2_s = new float[ WIDTH*HEIGHT*DEPTH ];
float *ute12_s = new float[ WIDTH*HEIGHT*DEPTH ];
float *ute12_air = new float[ WIDTH*HEIGHT*DEPTH ];

float *volume = new float[ WIDTH*HEIGHT*DEPTH ];

float *R2map = new float[ WIDTH*HEIGHT*DEPTH ];
float *R2map_low_LAC = new float[ WIDTH*HEIGHT*DEPTH ];
float *R2map_LAC = new float[ WIDTH*HEIGHT*DEPTH ];

float *R2bone_in_uteair = new float[ WIDTH*HEIGHT*DEPTH ];

float *umap_new = new float[ WIDTH*HEIGHT*DEPTH ];

int class_label = 0;
int *labels = new int[20000];
float *arr = new float[ WIDTH*HEIGHT*DEPTH ];
float *cluster_all = new float[ WIDTH*HEIGHT*DEPTH ];
float *cluster = new float[ WIDTH*HEIGHT*DEPTH ];

string atlasNII="";
string atlasMNC="";
string atlasmaskNII="";
string sinusAndNose="";
string basesinus="";
string baseNII="";
string sphenoid="";

string exec_with_output(const char* cmd) {
    char buffer[128];
    std::string result = "";
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    //cout << "exec output: " << result << endl;
    result.erase(std::remove(result.begin(), result.end(), '\n'), result.end());
    return result;
}

// Prepare binary command directories 
string str_dcm2mnc = exec_with_output("which dcm2mnc");
string str_nii2mnc = exec_with_output("which nii2mnc");
const char *dcm2mnc = str_dcm2mnc.c_str();
const char *nii2mnc = str_nii2mnc.c_str();

int _max( float const* a, size_t n )
{
	int maximum = 0;
	for( size_t i = 0; i < n; ++i )
		maximum = maximum < a[i] ? a[i] : maximum;
	return maximum;
}

float _mean( float const* a, size_t n )
{
	float m = 0.0f;
	for( size_t i = 0; i < n; ++i )
		m += a[i];
	return m/((float) n);
}

float find_delta( int* a, int* b)
{
	return sqrt( pow( a[0] - b[0] , 2) + pow( a[1] - b[1] , 2) );
}

float find_dist( int a, int b, int* c)
{
	return sqrt( pow( a - c[0] , 2) + pow( b - c[1] , 2) );
}

int* find_max( int* mean1, int* mean2, int m1, int m2)
{
	/* FIND LARGEST MEAN */
	int *largest = new int[2];
	largest[0] = largest[1] = 0;
	if(find_dist(0,0,mean1) > find_dist(0,0,mean2)){
		largest = mean1;
	} else {
		largest = mean2;
	}
	return largest;
}

int* find_mean(float* a,int m1,int m2,int label)
{
	int *new_mean = new int[2];
	new_mean[0] = 0;
	new_mean[1] = 0;
	int count = 0;
	for(size_t i = 0; i < m1; ++i)
	{
		for(size_t j = 0; j < m2; ++j)
		{
			if(a[ i + j*m1 ] == label){
				new_mean[0] += i;
				new_mean[1] += j;
				count += 1;
			}
		}
	}
	new_mean[0] = new_mean[0]/count;
	new_mean[1] = new_mean[1]/count;
	return new_mean;
}

void e_step(size_t m1, size_t m2, int* mean1, int* mean2)
{
	for(size_t i = 0; i < m1; ++i)
	{
		for(size_t j = 0; j < m2; ++j)
		{
			if(points[ i + j*m1 ] > 0){
				float dist1 = find_dist(i,j,mean1);
				float dist2 = find_dist(i,j,mean2);
				points[ i + j*m1 ] = (dist1 > dist2) ? 2 : 1;					
			}
		}
	}
}

inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int mni_register_brain_from_atlas(const char* dcmfolder){

	DIR *dir;
	struct dirent *ent;
	char *dcmfile;
	if((dir = opendir(dcmfolder)) != NULL){
		while((ent = readdir(dir)) != NULL){
			dcmfile = ent->d_name;
			if(dcmfile[0] != '.')
				break;
		}
		closedir(dir);
	}
	string dcmfile_full = string(dcmfolder) + string(dcmfile);

	cout << "Starting ANTS registration from MNI to patient space (about 2 minutes).." << endl;

	/* Get AGE of patient */
	int age;
	DcmFileFormat fileformat;
  	if (fileformat.loadFile(dcmfile_full.c_str()).good()){
  		const char *value = NULL;
  		DcmDataset *dataset = fileformat.getDataset();
  		if (dataset->findAndGetString(DCM_PatientAge, value, 4).good()){
  			string svalue(value);
  			string::size_type sz;
  			int age = stoi (svalue,&sz);
  			if (age <= 4){
  				cout << " - CHOSE ATLAS: 33-44w = 3y (should use more precise estimation of age!)" << endl;
  				atlasMNC="/usr/local/share/RESOLUTE/nihpd_asym_33-44_t1w.mnc";
  				atlasNII="/usr/local/share/RESOLUTE/nihpd_asym_33-44_t1w.nii";
				atlasmaskNII="/usr/local/share/RESOLUTE/nihpd_asym_33-44_mask.nii";
				sinusAndNose="/usr/local/share/RESOLUTE/nihpd_asym_33-44_sinus_and_nose2.nii";
				basesinus="/usr/local/share/RESOLUTE/nihpd_asym_33-44_basesinus.nii";
				baseNII="/usr/local/share/RESOLUTE/nihpd_asym_33-44_base11.nii";
				sphenoid = "/usr/local/share/RESOLUTE/nihpd_asym_33-44_nasal_septa_and_eithmoidal_sinuses.nii";
  			} else if(age < 8){
  				cout << " - CHOSE ATLAS: 4.5-8.5y" << endl;
  				atlasMNC="/usr/local/share/RESOLUTE/nihpd_asym_04.5-08.5_t1w.mnc";
  				atlasNII="/usr/local/share/RESOLUTE/nihpd_asym_04.5-08.5_t1w.nii";
				atlasmaskNII="/usr/local/share/RESOLUTE/nihpd_asym_04.5-08.5_mask.nii";
				sinusAndNose="/usr/local/share/RESOLUTE/nihpd_asym_04.5-08.5_sinus_and_nose2.nii";
				basesinus="/usr/local/share/RESOLUTE/nihpd_asym_04.5-08.5_basesinus.nii";
				baseNII="/usr/local/share/RESOLUTE/nihpd_asym_04.5-08.5_base11.nii";
				sphenoid = "/usr/local/share/RESOLUTE/nihpd_asym_04.5-08.5_nasal_septa_and_eithmoidal_sinuses.nii";
  			} else if(age <= 11){
  				cout << " - CHOSE ATLAS: 7-11y" << endl;
  				atlasMNC="/usr/local/share/RESOLUTE/nihpd_asym_07.0-11.0_t1w.mnc";
  				atlasNII="/usr/local/share/RESOLUTE/nihpd_asym_07.0-11.0_t1w.nii";
				atlasmaskNII="/usr/local/share/RESOLUTE/nihpd_asym_07.0-11.0_mask.nii";
				sinusAndNose="/usr/local/share/RESOLUTE/sinus_and_nose2.nii";
				basesinus="/usr/local/share/RESOLUTE/nihpd_asym_07.0-11.0_basesinus.nii";
				baseNII="/usr/local/share/RESOLUTE/nihpd_asym_07.0-11.0_base11.nii";
				sphenoid = "/usr/local/share/RESOLUTE/nihpd_asym_07.0-11.0_nasal_septa_and_eithmoidal_sinuses.nii";
  			} else if(age <= 18){
  				cout << " - CHOSE ATLAS: 13-18.5y" << endl;
  				atlasMNC="/usr/local/share/RESOLUTE/nihpd_asym_13.0-18.5_t1w.mnc";
  				atlasNII="/usr/local/share/RESOLUTE/nihpd_asym_13.0-18.5_t1w.nii";
				atlasmaskNII="/usr/local/share/RESOLUTE/nihpd_asym_13.0-18.5_mask.nii";
				sinusAndNose="/usr/local/share/RESOLUTE/nihpd_asym_13.0-18.5_sinus_and_nose2.nii";
				basesinus="/usr/local/share/RESOLUTE/nihpd_asym_13.0-18.5_basesinus.nii";
				baseNII="/usr/local/share/RESOLUTE/nihpd_asym_13.0-18.5_base11.nii";
				sphenoid = "/usr/local/share/RESOLUTE/nihpd_asym_13.0-18.5_nasal_septa_and_eithmoidal_sinuses.nii";
  			} else if (age > 18){
				cout << " - CHOSE ATLAS: >18y" << endl;
				atlasMNC="/usr/local/share/RESOLUTE/mni_icbm152_t1_tal_nlin_sym_09a.mnc";
				atlasNII="/usr/local/share/RESOLUTE/mni_icbm152_t1_tal_nlin_sym_09a.nii";
				atlasmaskNII="/usr/local/share/RESOLUTE/mni_icbm152_t1_tal_nlin_sym_09a_mask.nii";
				sinusAndNose="/usr/local/share/RESOLUTE/sinus_and_nose2.nii";
				basesinus="/usr/local/share/RESOLUTE/basesinus.nii";
				baseNII="/usr/local/share/RESOLUTE/base11.nii";
				sphenoid = "/usr/local/share/RESOLUTE/nasal_septa_and_eithmoidal_sinuses.nii";
			} else 
    			cout << "under 18!" << endl;
    	}
	} else {
		cout << "Failed to load DCM file.. stopping" << endl;
		return 1;
	}

	const string executable = string("ANTS 3 -m CC[") 
							+ string(atlasNII) 
							+ string(",/tmp/resolute_tmp/ute2.nii,1,4]")
							+ string(" -i 10x5x2 -o cc10x5x2 -t SyN[0.5] -r Gauss[3,0] -G")
							+ string(" >> /tmp/resolute_tmp/log_ants.txt");
	
	if(!file_exists("/tmp/resolute_tmp/cc10x5x2Affine.txt")){
		int status;
		status = system(executable.c_str());
		system("mv cc10x* /tmp/resolute_tmp/");
		cout << " - Finished ANTS.." << endl;
	} else {
		cout << " - ANTS already performed.." << endl; 
	}

}

int* find_scale_constants(){

	int max1 = _max(ute1,192*192*192);
	int max2 = _max(ute2,192*192*192);
	float *jhist = new float[ (max2+1) * (max1+1) ];
	for(size_t i = 0; i < max2 * max1; ++i)
	{
		jhist[i] = 0;
	}

	for(size_t i = 0; i < 192*192*192; i++)
	{
				int u_1 = ute1[ i ];
				int u_2 = ute2[ i ];
				jhist[ u_1 + u_2*max1 ]  += 1.0f;
	}
	float m = _mean(jhist,max1*max2);
	int *mean1 = new int[2];
	mean1[0] = 10;
	mean1[1] = 10;
	int *mean2 = new int[2];
	mean2[0] = (int) max2/4;
	mean2[1] = (int) max1/4;
	int *oldmean1 = new int[2];
	int *oldmean2 = new int[2];
	points = new float[ (max1+1) * (max2+1) ];
	for(size_t i = 0; i < max1; ++i)
	{
		for(size_t j = 0; j < max2; ++j)
		{
			points[ i + j*max1 ] = jhist[ i + j*max1 ] > 5*m ? 1 : 0;
		}
	}

	while(find_delta(mean1,oldmean1) > 5 || find_delta(mean2,oldmean2) > 5){
		// E STEP 
		e_step(max1,max2,mean1,mean2);	
		// M STEP 
		oldmean1 = mean1;
		oldmean2 = mean2;
		mean1 = find_mean(points,max1,max2,1);
		mean2 = find_mean(points,max1,max2,2);
	}
	int* largest = find_max(mean1,mean2,max1,max2);
	delete[] points;
	delete[] jhist;

	return largest;
}

void warp_image(const char* moving, const char* stationary, const char* new_name){
	const string executable = string("WarpImageMultiTransform 3 ") 
							+ string(moving) 
							+ string(" ")
							+ string(new_name)
							+ string(" -R ")
							+ string(stationary)
							+ string(" -i /tmp/resolute_tmp/cc10x5x2Affine.txt /tmp/resolute_tmp/cc10x5x2InverseWarp.nii.gz")
							+ string(" >> /tmp/resolute_tmp/log_warp.txt");
		
	if(file_exists("/tmp/resolute_tmp/cc10x5x2Affine.txt")){
		int status;
		status = system(executable.c_str());
	} else {
		cout << "Missing transformation file" << endl; 
	}
}

void system_call(vector<string> arg_input, const char *logfile){
	cout.flush();
	pid_t parent = getpid();
	pid_t pid = fork();
	if(pid > 0){
		int status;
		waitpid(pid, &status, 0);
	} else {		
		
		int fd; 
		if((fd = open(logfile, O_RDWR | O_CREAT))==-1){
  			perror("open");
		}
		dup2(fd,STDOUT_FILENO); 
		dup2(fd,STDERR_FILENO); 
		close(fd);

    	const char **argv = new const char* [arg_input.size()+1];
    	for (int j = 0;  j < arg_input.size();  j++){
            argv [j] = arg_input[j] .c_str();
    	}
    	argv [arg_input.size()] = NULL;

		execv(arg_input[0].c_str(),(char **)argv);
	}
}

int* get_location(int i){
	int *out = new int[3];
	int slice = floor( i / (192*192) );
	int row = floor( ( i - 192*192*slice ) / 192 );
	int col = i - 192*192*slice - 192*row;
	out[0] = slice;
	out[1] = row;
	out[2] = col;
	return out;
}

float blur_voxel(int i, float *arr){

	int *loc_i = get_location(i);

	int i_tl = i-192-1;
	int i_tr = i-192+1;
	int i_t = i-192;
	int i_l = i-1;
	int i_r = i+1;
	int i_bl = i+192-1;
	int i_b = i+192;
	int i_br = i+192+1;
	int *loc_i_tl = get_location(i_tl);
	int *loc_i_br = get_location(i_br);
	if(loc_i[0] == loc_i_tl[0] && loc_i[0] == loc_i_br[0] // same slice
		&& loc_i[1] == (loc_i_tl[1]+1) && loc_i[1] == (loc_i_br[1]-1) // 1 row diff
		&& loc_i[2] == (loc_i_tl[2]+1) && loc_i[2] == (loc_i_br[2]-1) // 1 col diff
		&& loc_i[1] > 0 && loc_i[2] > 0) // top left corners
	{
		return 0.01*arr[i_tl]+0.08*arr[i_t]+0.01*arr[i_tr]
			 + 0.08*arr[i_l]+0.64*arr[i]+0.08*arr[i_r]
		  	 + 0.01*arr[i_bl]+0.08*arr[i_b]+0.01*arr[i_br];
	} else {
		return arr[i];
	}
}

void load_raw_files(){
	ifstream in( "/tmp/resolute_tmp/ute1.raw", ios::in | ios::binary );
	in.read( reinterpret_cast< char* >( ute1 ), sizeof(float)*WIDTH*HEIGHT*DEPTH );
	in.close();
	ifstream in2( "/tmp/resolute_tmp/ute2.raw", ios::in | ios::binary );
	in2.read( reinterpret_cast< char* >( ute2 ), sizeof(float)*WIDTH*HEIGHT*DEPTH );
	in2.close();

	// Warp the sinus and nose area
	warp_image(sinusAndNose.c_str(), "/tmp/resolute_tmp/ute2.nii", "/tmp/resolute_tmp/sinus.nii");
	vector<string> vector_sinusAndNose {nii2mnc,"/tmp/resolute_tmp/sinus.nii","/tmp/resolute_tmp/sinus2.mnc"};
	system_call(vector_sinusAndNose,"/tmp/resolute_tmp/log_trash.txt");
	system("mincresample /tmp/resolute_tmp/sinus2.mnc -like /tmp/resolute_tmp/ute2.mnc /tmp/resolute_tmp/sinusmask.mnc -clobber -quiet");
	system("minctoraw -nonormalize -float /tmp/resolute_tmp/sinusmask.mnc > /tmp/resolute_tmp/sinusmask.raw");
	ifstream in_sinus( "/tmp/resolute_tmp/sinusmask.raw", ios::in | ios::binary );
	in_sinus.read( reinterpret_cast< char* >( sinusmask ), sizeof(float)*WIDTH*HEIGHT*DEPTH );
	in_sinus.close();

	// Extract brain and CSF values
	cout << "Extracting brain and CSF values (about 30 seconds)" << endl;
	warp_image(atlasmaskNII.c_str(), "/tmp/resolute_tmp/ute2.nii", "/tmp/resolute_tmp/atlasbrain.nii");
	vector<string> vector_atlasbrain {nii2mnc,"/tmp/resolute_tmp/atlasbrain.nii","/tmp/resolute_tmp/atlasbrain.mnc"};
	system_call(vector_atlasbrain,"/tmp/resolute_tmp/log_trash.txt");
	system("mincresample /tmp/resolute_tmp/atlasbrain.mnc -like /tmp/resolute_tmp/ute2.mnc /tmp/resolute_tmp/brainmask.mnc -quiet -clobber");
	system("mincmath -quiet -gt -const 0 /tmp/resolute_tmp/brainmask.mnc /tmp/resolute_tmp/brain.mnc -clobber");
	system("minctoraw -nonormalize -float /tmp/resolute_tmp/brain.mnc > /tmp/resolute_tmp/brain.raw");
	ifstream in_brain( "/tmp/resolute_tmp/brain.raw", ios::in | ios::binary );
	in_brain.read( reinterpret_cast< char* >( brain ), sizeof(float)*WIDTH*HEIGHT*DEPTH );
	in_brain.close();
	// CSF 
	if(!file_exists("/tmp/resolute_tmp/ute2_nuc.mnc")){
		vector<string> vector_insect {"/usr/local/share/RESOLUTE/insect_modified","/tmp/resolute_tmp/ute2.mnc",atlasMNC};
		system_call(vector_insect,"/tmp/resolute_tmp/log_insect.txt");
		system("mv ute2_* /tmp/resolute_tmp/");
	}
	system("mincresample /tmp/resolute_tmp/ute2_nuc_tal_cla.mnc -like /tmp/resolute_tmp/ute2.mnc -transformation /tmp/resolute_tmp/ute2_nuc_total.xfm -invert_transformation /tmp/resolute_tmp/ute2_cla.mnc -quiet -clobber");
	system("mincresample -quiet /tmp/resolute_tmp/ute2_nuc_tal_cla.mnc -invert_transformation -transformation /tmp/resolute_tmp/ute2_nuc_total.xfm /tmp/resolute_tmp/ute2_nuc_tal_cla_rsl.mnc -like /tmp/resolute_tmp/ute2.mnc -clobber");
	system("mincmath -quiet -eq -const 1 /tmp/resolute_tmp/ute2_nuc_tal_cla_rsl.mnc /tmp/resolute_tmp/csf.mnc -clobber");
	system("minctoraw -nonormalize -float /tmp/resolute_tmp/csf.mnc > /tmp/resolute_tmp/csf.raw");
	ifstream in_csf( "/tmp/resolute_tmp/csf.raw", ios::in | ios::binary );
	in_csf.read( reinterpret_cast< char* >( csf ), sizeof(float)*WIDTH*HEIGHT*DEPTH );
	in_csf.close();
	cout << " - Done" << endl;

	// Extract base of skull mask
	cout << "Running extract masks" << endl;
	warp_image(baseNII.c_str(),"/tmp/resolute_tmp/ute2.nii","/tmp/resolute_tmp/base.nii");
	vector<string> vector_base {nii2mnc,"/tmp/resolute_tmp/base.nii","/tmp/resolute_tmp/base.mnc"};
	system_call(vector_base,"/tmp/resolute_tmp/log_trash.txt");
	system("mincresample /tmp/resolute_tmp/base.mnc -like /tmp/resolute_tmp/ute2.mnc /tmp/resolute_tmp/base2.mnc -clobber -quiet");
	system("mincmath -quiet -segment -const2 0.5 1.5 /tmp/resolute_tmp/base2.mnc /tmp/resolute_tmp/base.mnc -clobber");
	system("minctoraw -nonormalize -float /tmp/resolute_tmp/base.mnc > /tmp/resolute_tmp/base.raw");
	ifstream in_base( "/tmp/resolute_tmp/base.raw", ios::in | ios::binary );
	in_base.read( reinterpret_cast< char* >( base ), sizeof(float)*WIDTH*HEIGHT*DEPTH );
	in_base.close();

	warp_image(basesinus.c_str(), "/tmp/resolute_tmp/ute2.nii", "/tmp/resolute_tmp/basesinus.nii");
	vector<string> vector_basesinus {nii2mnc,"/tmp/resolute_tmp/basesinus.nii","/tmp/resolute_tmp/warped_basesinus.mnc"};
	system_call(vector_basesinus,"/tmp/resolute_tmp/log_trash.txt");
	system("mincresample /tmp/resolute_tmp/warped_basesinus.mnc -like /tmp/resolute_tmp/ute2.mnc /tmp/resolute_tmp/r2noise_rsl.mnc -clobber -quiet");
	system("mincmath -quiet -segment -const2 0.5 1.5 /tmp/resolute_tmp/r2noise_rsl.mnc /tmp/resolute_tmp/r2noise.mnc -clobber");
	system("minctoraw -nonormalize -float /tmp/resolute_tmp/r2noise.mnc > /tmp/resolute_tmp/r2noise.raw");
	ifstream in_r2noise( "/tmp/resolute_tmp/r2noise.raw", ios::in | ios::binary );
	in_r2noise.read( reinterpret_cast< char* >( r2noise ), sizeof(float)*WIDTH*HEIGHT*DEPTH );
	in_r2noise.close();

	warp_image(sphenoid.c_str(), "/tmp/resolute_tmp/ute2.nii", "/tmp/resolute_tmp/sphenoid.nii");
	vector<string> vector_sphenoid {nii2mnc,"/tmp/resolute_tmp/sphenoid.nii","/tmp/resolute_tmp/sphenoid.mnc"};
	system_call(vector_sphenoid,"/tmp/resolute_tmp/log_trash.txt");
	system("mincresample /tmp/resolute_tmp/sphenoid.mnc -like /tmp/resolute_tmp/ute2.mnc /tmp/resolute_tmp/sphenoid2.mnc -clobber -quiet");
	system("mincmath -quiet -gt -const 0 /tmp/resolute_tmp/sphenoid2.mnc /tmp/resolute_tmp/sphenoid.mnc -clobber");
	system("minctoraw -nonormalize -float /tmp/resolute_tmp/sphenoid.mnc > /tmp/resolute_tmp/sphenoid.raw");
	ifstream in_sphenoid( "/tmp/resolute_tmp/sphenoid.raw", ios::in | ios::binary );
	in_sphenoid.read( reinterpret_cast< char* >( sphenoidmask ), sizeof(float)*WIDTH*HEIGHT*DEPTH );
	in_sphenoid.close();
	cout << " - Done" << endl;
}

void scale_utes(){
	int *scales = find_scale_constants();
	//cout <<  scales[0] << " " << scales[1] << endl;

	for(size_t i=0; i<WIDTH*HEIGHT*DEPTH; i++){
		ute1_s[i] = ute1[i]/scales[0]*1000;
		ute2_s[i] = ute2[i]/scales[1]*1000;
		ute12_s[i] = ute1_s[i]+ute2_s[i];

		R2map[i] = ( log(ute1[i])-log(ute2[i]) ) / 2.39;
	}
}

void update_labels(int from, int to){
	for(size_t l=0; l<class_label; l++){
		if(labels[l] == from)
			labels[l] = to;
	}
}

int cluster_volume(int i){

	int *loc_i = get_location(i);

	int i_t = i-WIDTH;
	int i_l = i-1;
	int i_p = i-WIDTH*HEIGHT;

	int *loc_l = get_location(i_l);
	int *loc_t = get_location(i_t);
	int *loc_p = get_location(i_p);

	if(i==0){
		class_label++;
		labels[class_label] = class_label;
		return class_label;
	} else {

		bool update = 0;

		// Not air
		if(arr[i] > 0.5)
			return 0;

		// Same as left neighbour
		else if(loc_i[2]>0 && arr[i] == arr[i_l]){ // not on border (no left)

			// Check other neighbours
			if(loc_i[1]>0 && arr[i] == arr[i_t]){

				update = 1;

				// Same as top, update pointer to point to lowest number
				if(labels[(int)cluster_all[i_t]] > labels[(int)cluster_all[i_l]] && labels[(int)cluster_all[i_l]] > 0){
					int old_label = labels[(int)cluster_all[i_t]];
					labels[(int)cluster_all[i_t]] = (labels[(int)cluster_all[i_l]]);
					update_labels(old_label,labels[(int)cluster_all[i_l]]);
				}
				if(labels[(int)cluster_all[i_t]] < labels[(int)cluster_all[i_l]] && labels[(int)cluster_all[i_t]] > 0){
					int old_label = labels[(int)cluster_all[i_l]];
					labels[(int)cluster_all[i_l]] = (labels[(int)cluster_all[i_t]]);
					update_labels(old_label,labels[(int)cluster_all[i_t]]);
				}

			} 

			if(loc_i[0] > 0 && arr[i] == arr[i_p]){

				// Same as previous slice, update pointer to point to lowest number
				if(labels[(int)cluster_all[i_p]] > labels[(int)cluster_all[i_l]] && labels[(int)cluster_all[i_l]] > 0){
					int old_label = labels[(int)cluster_all[i_p]];
					labels[(int)cluster_all[i_p]] = (labels[(int)cluster_all[i_l]]);
					update_labels(old_label,labels[(int)cluster_all[i_l]]);
				}
				if(labels[(int)cluster_all[i_p]] < labels[(int)cluster_all[i_l]] && labels[(int)cluster_all[i_p]] > 0){
					int old_label = labels[(int)cluster_all[i_l]];
					labels[(int)cluster_all[i_l]] = (labels[(int)cluster_all[i_p]]);
					if(update)
						update_labels(labels[(int)cluster_all[i_t]],labels[(int)cluster_all[i_p]]);
					update_labels(old_label,labels[(int)cluster_all[i_p]]);
				}

			}
			
			return labels[(int)cluster_all[i_l]];
		}

		// Same as top neighbour
		else if(loc_i[1]>0 && arr[i] == arr[i_t]){ // not on first row
			
			// Check other neighbours
			// Dont need to check left, otherwise it would have been updated in last if.
			if(loc_i[0] > 0 && arr[i] == arr[i_p]){
				// Same as previous slice, update pointer to point to lowest number
				if(labels[(int)cluster_all[i_p]] > labels[(int)cluster_all[i_t]] && labels[(int)cluster_all[i_t]] > 0){
					int old_label = labels[(int)cluster_all[i_p]];
					labels[(int)cluster_all[i_p]] = (labels[(int)cluster_all[i_t]]);
					update_labels(old_label,labels[(int)cluster_all[i_t]]);
				}
				if(labels[(int)cluster_all[i_p]] < labels[(int)cluster_all[i_t]] && labels[(int)cluster_all[i_p]] > 0){
					int old_label = labels[(int)cluster_all[i_t]];
					labels[(int)cluster_all[i_t]] = (labels[(int)cluster_all[i_p]]);
					update_labels(old_label,labels[(int)cluster_all[i_p]]);
				}
			}

			return labels[(int)cluster_all[i_t]];
		}

		// Same as previous neighbour
		else if(loc_i[0] > 0 && arr[i] == arr[i_p]){
			// No need to check other neighbours, top and left have different values.
			return labels[(int)cluster_all[i_p]];
		} 

		// No neighbours matched. We create a new color.
		else {
			class_label++;
			labels[class_label] = class_label;
			return class_label;
		}
	}
}

void locate_inner_air(){
	cout << "Calculating inner air using cluster location" << endl;
	for(size_t i = 0; i < WIDTH*HEIGHT*DEPTH; i++){
		arr[i] = (blur_voxel(i,ute12_s) < 600) ? 0.0 : 1.0;
		cluster_all[i] = cluster_volume(i);
	}

	int actual_clusters = 0;
	int *clusters_arr = new int[class_label];
	for(size_t i = 0; i < class_label; i++)
		clusters_arr[i] = 0;
	for(size_t i = 0; i < WIDTH*HEIGHT*DEPTH; i++){
		cluster_all[i] = labels[(int)cluster_all[i]];

		if( (int)cluster_all[i] == 0)
			continue;
		else if( clusters_arr[ (int)cluster_all[i] ] == 0 ){
			actual_clusters++;
			clusters_arr[ (int)cluster_all[i] ] = actual_clusters;
			cluster[i] = actual_clusters;
		} else {
			cluster[i] = clusters_arr[ (int)cluster_all[i] ];
		}
	}
	delete[] clusters_arr;
	delete[] arr;
	delete[] labels;

	cout << " - Found " << actual_clusters << " clusters." << endl;
}

void calculate_umap(){
	cout << "Calculating tissue classes and combining masks" << endl;
	for(size_t i=0; i<WIDTH*HEIGHT*DEPTH; i++){

		/* Calculate patient volume and air regions */
		if(blur_voxel(i,ute12_s) < 600){
			ute12_air[i] = 1.0;
			volume[i] = 0.0;
		} else {
			volume[i] = 1.0;
			ute12_air[i] = 0.0;
		}

		/* Calculate R2star maps */
		float vox_R2blur = blur_voxel(i,R2map);
		float vox_R2s_max = (vox_R2blur > R2map[i]) ? vox_R2blur*1000.0 : R2map[i]*1000.0;
		float vox_HU = 0.000001351*pow(vox_R2s_max,3.0) - 0.003617*pow(vox_R2s_max,2.0) + 3.841*vox_R2s_max - 19.46;
		float vox_LAC_tmp = (0.000051*(vox_HU+1000)+0.0471)*10000;
		float vox_LAC = (vox_R2s_max > 100 && vox_LAC_tmp > 1010) ? vox_LAC_tmp : 0.0;
		R2map_LAC[i] = vox_LAC;

		float vox_low_HU = 0.000001351*pow(R2map[i],3.0) - 0.003617*pow(R2map[i],2.0) + 3.841*R2map[i] - 19.46;
		float vox_low_LAC_tmp = (0.000051*(vox_low_HU+1000)+0.0471)*10000;
		float vox_low_LAC = (R2map[i] > 300 && vox_low_LAC_tmp > 1010) ? vox_low_LAC_tmp : 0.0;
		R2map_low_LAC[i] = vox_low_LAC;

		/* Extract ute signal in air */
		if(vox_R2blur > 0.3 && ute12_air[i] > 0.5 && sinusmask[i] < 0.5 && cluster[i] > 1.5)
			R2bone_in_uteair[i] = 1.0;
		else
			R2bone_in_uteair[i] = 0.0;

		
		/* Combine the info */
        if(volume[i] > 0.5){
        	// Inside volume
        	if(brain[i] > 0.5 && R2map_LAC[i] < 1500){
        		if(csf[i] > 0.5)
        			umap_new[i] = 960;
				else
					umap_new[i] = 990;
         	} else if(R2map_LAC[i] > 1010){
         		if(r2noise[i] > 0.5)
         			umap_new[i] = 1100;
         		else if(base[i] > 0.5 && R2map_low_LAC[i] > 1010)
         			umap_new[i] = R2map_low_LAC[i];
         		else if(base[i] > 0.5 || ute2_s[i] > 1200)
         			umap_new[i] = 925;
         		else
         		    umap_new[i] = R2map_LAC[i];
        	} else {

        		/* find nesal septa and ethmoidal sinuses */
        		if(ute12_s[i] < 800 && sphenoidmask[i] > 0.5)
        			umap_new[i] = 100;
        		else if (ute12_s[i] < 1600 && sphenoidmask[i] > 0.5)
        			umap_new[i] = 600;
        		else	
	        		umap_new[i] = 940;
        	}
        } else if(R2bone_in_uteair[i] > 0.5){
        	umap_new[i] = 1000;
        } else {
        	umap_new[i] = 0;
        }

	}

}

void prepare_mnc_and_nifty_files(const char *argv[]){
	cout << "Building .mnc files" << endl;

	system("mkdir /tmp/resolute_tmp");
	system("touch /tmp/resolute_tmp/log_trash.txt");

	vector<string> vector_ute2 {dcm2mnc,argv[2],"-dname","","-fname","ute2","/tmp/resolute_tmp/","-clobber"};
	system_call(vector_ute2,"/tmp/resolute_tmp/log_trash.txt");
	cout << " - ute2 done" << endl;

	vector<string> vector_ute1 {dcm2mnc,argv[1],"-dname","","-fname","ute1","/tmp/resolute_tmp/","-clobber"};
	system_call(vector_ute1,"/tmp/resolute_tmp/log_trash.txt");
	cout << " - ute1 done" << endl;

	cout << "Building .nii files" << endl;
	system("mnc2nii /tmp/resolute_tmp/ute2.mnc /tmp/resolute_tmp/ute2.nii");

	cout << "Building .raw files" << endl;
	system("minctoraw -nonormalize -float -unsigned /tmp/resolute_tmp/ute1.mnc > /tmp/resolute_tmp/ute1.raw");
	system("minctoraw -nonormalize -float -unsigned /tmp/resolute_tmp/ute2.mnc > /tmp/resolute_tmp/ute2.raw");
}

void save_to_dcm(const char* uteumapfolder, const char *out_folder){
	string executable_mkdir = "mkdir "+string(out_folder);
	system(executable_mkdir.c_str());

	DIR *dir;
	struct dirent *ent;
	char *dcmfile;
	if((dir = opendir(uteumapfolder)) != NULL){
		while((ent = readdir(dir)) != NULL){
			dcmfile = ent->d_name;
			if(strlen(dcmfile) >= 4 && ( strcmp(dcmfile+strlen(dcmfile)-4,".dcm") == 0
									  || strcmp(dcmfile+strlen(dcmfile)-4,".DCM") == 0
									  || strcmp(dcmfile+strlen(dcmfile)-4,".IMA") == 0
									  || strcmp(dcmfile+strlen(dcmfile)-4,".ima") == 0))
			{

				string dcmfile_full = string(uteumapfolder) + string(dcmfile);
				string filetype = string(dcmfile+strlen(dcmfile)-3);
				int instance_number;
				DcmFileFormat fileformat;
  				if (fileformat.loadFile(dcmfile_full.c_str()).good()){
  					int value = 0;
  					DcmDataset *dataset = fileformat.getDataset();
  					if (dataset->findAndGetSint32(DCM_InstanceNumber, value).good()){

  						char executable[200];
						snprintf(executable, sizeof(executable), "cp %s %s/IM-0001-%04d.%s",dcmfile_full.c_str(),out_folder,static_cast<int>(value),filetype.c_str());
						system(executable);

					}
				}
			}
		}
		closedir(dir);
	}

    for(size_t d = 0; d < DEPTH; d++){

    	Uint16 *slice = new Uint16[WIDTH*HEIGHT];
    	for(size_t j = 0; j < HEIGHT; j++){
        	for(size_t i = 0; i < WIDTH; i++){
				slice[i + j*HEIGHT] = umap_new[i + j*HEIGHT + d*WIDTH*HEIGHT];
	      	}
	    }
	    DcmFileFormat fileformat;
		char buff2[100];
		snprintf(buff2, sizeof(buff2), "%s/IM-0001-%04d.dcm",out_folder,static_cast<int>(d+1));
		const char* name2 = buff2;	      			
		if (fileformat.loadFile(name2).good()){
	        DcmDataset *dataset = fileformat.getDataset();
	        dataset->putAndInsertUint16Array(DCM_PixelData, slice, 192*192);
	        dataset->putAndInsertString(DCM_SeriesDescription, "RESOLUTE");
	        dataset->putAndInsertUint16(DCM_SeriesNumber, 100);

	        OFCondition status = fileformat.saveFile(name2, EXS_LittleEndianExplicit);
	        if(status.bad())
	        	cerr << "Error: " << status.text() << endl;
		}
		delete[] slice;
	}

	cout << "Finished with patient. Saved data to: " << out_folder << endl;
}

int main(int argc, const char *argv[]) {

	if(argc < 5){
		cout << "Error running RESOLUTE. Must supply exactly 4 inputs, you supplied " << argc << endl;
		return -1;
	}

	cout << "Running RESOLUTE with parameters " << endl;
	cout << "\tUTE TE1 folder: " << argv[1] << endl;
	cout << "\tUTE TE2 folder: " << argv[2] << endl;
	cout << "\tUTE Umap folder: " << argv[3] << endl;
	cout << "\tOutput DCM dir: " << argv[4] << endl; 
	
	/* Create .mnc, .nii and .dcm files */
	prepare_mnc_and_nifty_files(argv);

	/* Align to ICBM atlas */
	mni_register_brain_from_atlas(argv[1]);

	/* Load raw files */
	load_raw_files();

	/* Scale UTEs */
	scale_utes();

	/* Find air inside patient */
	locate_inner_air();

	/* Calculate tissue maps and combine into umap */
	calculate_umap();

	/* Save data to DCM */
	save_to_dcm(argv[3], argv[4]);

	system("rm -rf /tmp/resolute_tmp");
	
	delete[] ute1;
	delete[] ute2;
	delete[] sinusmask;
	delete[] brain;
	delete[] csf;
	delete[] base;
	delete[] r2noise;
	delete[] sphenoidmask;
	delete[] ute1_s;
	delete[] ute2_s;
	delete[] ute12_s;
	delete[] ute12_air;
	delete[] volume;
	delete[] R2map;
	delete[] R2map_low_LAC;
	delete[] R2map_LAC;
	delete[] R2bone_in_uteair;
	delete[] umap_new;
	delete[] cluster;

	return 0;

}