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

// Pour prendre une tranche d'une matrice existante qui partage la même mémoire.
float** matrix_view(float** mat, int i, int j, int length, int width) {
    float ** view = malloc(sizeof(*view) * length);

    for (int v = 0;v < length;v++) {
        view[v] = &mat[v + i][j];
    }

    return view;
}

float ** denoise(float T, float ** mat, int length, int width) {
    float ** ImgDenoised_temp = fmatrix_allocate_2d(length, width);
    fmatrix_move(length, width, mat, ImgDenoised_temp);

    float ** images[8*8];

    for (int decali = 0; decali < 8;decali++) {
        for (int decalj = 0; decalj < 8;decalj++) {
            decal_toroid(ImgDenoised_temp, length, width, decali, decalj);
            for (int blocki = 0;8*blocki < length;blocki++) {
                for (int blockj = 0;8*blockj < width;blockj++) {
                    float ** view = matrix_view(ImgDenoised_temp, 8*blocki,  8*blockj, 8, 8);
                    ddct8x8s(-1, view);
                    for (int i = 0;i < 8;i++) {
                        for (int j = 0;j < 8;j++) {
                            float val = view[i][j]; 
                            if (val*val < T*T)
                                view[i][j] = 0.0;
                        }                     
                    }
                    ddct8x8s(1, view);
                    free(view);
                }             
            }
            decal_toroid(ImgDenoised_temp, length, width, length-decali, width-decalj);
            int ind = decali + decalj*8;
            // printf("i: %d, j: %d, ind: %d\n", decali, decalj, ind);
            images[ind] = fmatrix_allocate_2d(length, width);
            fmatrix_move(length, width, ImgDenoised_temp, images[decali + decalj*8]);
            printf("decal[%d,%d] done\n", decali, decalj);
        }
    }

    float ** moyenne = fmatrix_allocate_2d(length, width);
    fmatrix_zero(length, width, moyenne);
    for (int i = 0;i < length;i++) {
        for (int j = 0;j < width;j++) {
            for (int n = 0;n < 8*8;n++) {
                moyenne[i][j] += images[n][i][j];
            }
            moyenne[i][j] /= 64;
        }
    }
    for (int i = 0;i < 8*8;i++) {
        free_fmatrix_2d(images[i]);
    }
    free_fmatrix_2d(ImgDenoised_temp);
    return moyenne;
}



//---------------------------------------------------------
//---------------------------------------------------------
// PROGRAMME PRINCIPAL   ----------------------------------                     
//---------------------------------------------------------
//---------------------------------------------------------
int main(int argc,char** argv)
{
 int length,width;
 int bestlim =0;
 char BufSystVisuImg[NBCHAR];

 //>Lecture Image 
 float** Img=LoadImagePgm(NAME_IMG_IN,&length,&width);
 
 //>Allocation memory
 float** ImgDegraded=fmatrix_allocate_2d(length,width);
 float** ImgDenoised=fmatrix_allocate_2d(length,width);  


 //>Degradation 
 copy_matrix(ImgDegraded,Img,length,width);
 add_gaussian_noise(ImgDegraded,length,width,SIGMA_NOISE*SIGMA_NOISE);

 printf("Meilleure Limite: %.i\n",bestlim);
 printf("\n\n  Info Noise");
 printf("\n  ------------");
 printf("\n    > MSE = [%.2f]",computeMMSE(ImgDegraded,Img,length)); 
 

 //=========================================================
 //== PROG =================================================
 //=========================================================



 float best_T = -1;
 float best_error = 1e10;
 float ** best_image;

 for (float T = 10.0;T < 25.0;T += 1.0) {
     float ** denoised = denoise(T, ImgDegraded, length, width);
     float err = computeMMSE(Img, denoised, length);
     if (err < best_error) {
         best_error = err;
         best_T = T;
         best_image = denoised;
     }
 }

 printf("best T: %f\n", best_T);
 printf("best error: %f\n", best_error);
  //---------------------------------------------
  // SAUVEGARDE 
  // -------------------
  // L'image d�grad�e             > ImgDegraded
  // Le resultat du debruitage    > ImgFiltered
  //----------------------------------------------
  SaveImagePgm(NAME_IMG_DEG,ImgDegraded,length,width);
  SaveImagePgm(NAME_IMG_OUT,best_image,length,width);  

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
