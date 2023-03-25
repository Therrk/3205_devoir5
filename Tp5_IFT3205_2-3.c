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
  int i,j,k,l;
  float bestMMSE=999999999;
  int bestlim =0;
  int lim =0;
  float MMSE, total;
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

  do {
	  total=0;
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
  	        if (fabs(Imagette[k][l])<lim) {
              Imagette[k][l]=0;
            }
            total+=fabs(Imagette[k][l]);
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
  	  bestlim = lim;
  	  fmatrix_move(length, width, ImgDenoised_temp, ImgDenoised);
    }
    lim+=1;
  } while (total!=0);

 printf("Meilleure Limite: %.i\n",bestlim);
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
