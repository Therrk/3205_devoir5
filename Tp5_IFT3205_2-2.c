//------------------------------------------------------
// Prog    : Tp5_IFT3205                          
// Auteur  : Élie Leblanc, Justin veilleux
// Date    :                                  
// version :                                             
// langage : C                                          
// labo    : DIRO                                       
//------------------------------------------------------
// elie.leblanc@umontreal.ca justin.veilleux@umontreal.ca


//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include "FonctionDemo5.h"
#include "fonctions.c"

//------------------------------------------------
// DEFINITIONS -----------------------------------
//------------------------------------------------
#define NAME_VISUALISER "display "
#define NAME_IMG_IN  "lena512"
#define NAME_IMG_OUT  "lena512_Restored"
#define NAME_IMG_DEG  "lena512_Degraded"

//------------------------------------------------
// PROTOTYPE DE FONCTIONS  -----------------------
//------------------------------------------------
//----------------------------------------------------------
// copy matrix
//----------------------------------------------------------
void copy_matrix(float** mat1,float** mat2,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) mat1[i][j]=mat2[i][j];
}



//------------------------------------------------
// CONSTANTES ------------------------------------
//------------------------------------------------
#define SIGMA_NOISE  30 

int coords[64][2] =  {{0,0},{0,1},{1,0},{2,0},{1,1},{0,2},{0,3},{1,2},{2,1},{3,0},{4,0},{3,1},{2,2},{1,3},{0,4},{0,5},{1,4},{2,3},{3,2},{4,1},{5,0},{6,0},{5,1},{4,2},{3,3},{2,4},{1,5},{0,6},{0,7},{1,6},{2,5},{3,4},{4,3},{5,2},{6,1},{7,0},{7,1},{6,2},{5,3},{4,4},{3,5},{2,6},{1,7},{2,7},{3,6},{4,5},{5,4},{6,3},{7,2},{7,3},{6,4},{5,5},{4,6},{3,7},{4,7},{5,6},{6,5},{7,4},{7,5},{6,6},{5,7},{7,6},{6,7},{7,7}};
//------------------------------------------------
//------------------------------------------------
// FONCTIONS  ------------------------------------   
//------------------------------------------------
//------------------------------------------------

//---------------------------------------------------------
//---------------------------------------------------------
// PROGRAMME PRINCIPAL   ----------------------------------                     
//---------------------------------------------------------
//---------------------------------------------------------
int main(int argc,char** argv)
{
 int length,width;
  int i,j,k,l,m;
  int bestM=0;
  float bestMMSE=999999999;
  float MMSE;
 char BufSystVisuImg[NBCHAR];

 //>Lecture Image 
 float** Img=LoadImagePgm(NAME_IMG_IN,&length,&width);
 
 //>Allocation memory
 float** ImgDegraded=fmatrix_allocate_2d(length,width);
 float** ImgDenoised=fmatrix_allocate_2d(length,width);  
 float** ImgDenoised_temp=fmatrix_allocate_2d(length,width);  
  float** Imagette = fmatrix_allocate_2d(8,8);

 //>Degradation 
 copy_matrix(ImgDegraded,Img,length,width);
 add_gaussian_noise(ImgDegraded,length,width,SIGMA_NOISE*SIGMA_NOISE);

for (m = 0; m < 64; m++) {
 for (i = 0; i < length/8; i++) {
	  for (j = 0; j < width/8; j++) {
	    for (k = 0; k < 8; k++) {
	      for (l = 0; l < 8; l++) {
	        Imagette[k][l]=ImgDegraded[(8*i)+k][(8*j)+l];
        }
      }
      ddct8x8s(1,Imagette);
      for (k = 0; k < 8; k++) {
	      for (l = 0; l < 8; l++) {
	        if ((8*k)+l>m) {
	          Imagette[coords[(8*k)+l][0]][coords[(8*k)+l][1]]=0;
          }
        }
      }
      ddct8x8s(-1,Imagette);
      for (k = 0; k < 8; k++) {
	      for (l = 0; l < 8; l++) {
	        ImgDenoised_temp[(8*i)+k][(8*j)+l]=Imagette[k][l];
        }
      }
    }
  }
  MMSE = computeMMSE(Img,ImgDenoised, length);
  if (MMSE<bestMMSE) {
	  bestMMSE = MMSE;
	  bestM = m;
	  fmatrix_move(length, width, ImgDenoised_temp, ImgDenoised);
  }
}
 
 printf("\n\n  Info Noise");
 printf("\n  ------------");
 printf("\n    > MSE = [%.2f]",computeMMSE(ImgDegraded,Img,length)); 
 

 //=========================================================
 //== PROG =================================================
 //=========================================================

 // .....

  //---------------------------------------------
  // SAUVEGARDE 
  // -------------------
  // L'image d�grad�e             > ImgDegraded
  // Le resultat du debruitage    > ImgFiltered
  //----------------------------------------------
  SaveImagePgm(NAME_IMG_DEG,ImgDegraded,length,width);
  SaveImagePgm(NAME_IMG_OUT,ImgDenoised,length,width);  

  //>Visu Img
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_IN);
  strcat(BufSystVisuImg,".pgm&");
  printf("\n > %s",BufSystVisuImg);
  system(BufSystVisuImg);

  //Visu ImgDegraded
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_DEG);
  strcat(BufSystVisuImg,".pgm&");
  printf("\n > %s",BufSystVisuImg);
  system(BufSystVisuImg);
  
  //Visu ImgFiltered
  strcpy(BufSystVisuImg,NAME_VISUALISER);
  strcat(BufSystVisuImg,NAME_IMG_OUT);
  strcat(BufSystVisuImg,".pgm&");
  printf("\n > %s",BufSystVisuImg);
  system(BufSystVisuImg);

    
//--------------- End -------------------------------------     
//----------------------------------------------------------

  //Liberation memoire pour les matrices
  if (Img)          free_fmatrix_2d(Img);
  if (ImgDegraded)  free_fmatrix_2d(ImgDegraded);
  if (ImgDenoised)  free_fmatrix_2d(ImgDenoised);
  
  //Return
  printf("\n C'est fini... \n");; 
  return 0;
 }
 


//----------------------------------------------------------
// Allocation matrix 3d float
// --------------------------
//
// float*** fmatrix_allocate_3d(int dsize,int vsize,int hsize)
// > Alloue dynamiquement un tableau 3D [dsize][vsize][hsize]
//
//-----------------------------------------------------------

//----------------------------------------------------------
//  ddct8x8s(int isgn, float** tab)
// ---------------------------------
//
// Fait la TCD (sgn=-1) sur un tableau 2D tab[8][8]
// Fait la TCD inverse (sgn=1) sur un tableau 2D tab[8][8]
//
//-----------------------------------------------------------
