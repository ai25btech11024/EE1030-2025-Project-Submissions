#include <stdio.h>
#include "matrix.h"
#include <stdlib.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main() {
	int k;
	printf("Enter the value of k: ");
	scanf("%d",&k);

	char of[100];
	char nf[100];

	printf("Enter name of input file: ");
	scanf("%s",of);

	printf("Enter name of output file: ");
	scanf("%s",nf);

	int w,h,max;
	unsigned char *img = stbi_load(of,&w,&h,&max,1);

	if(!img) {
		printf("Could not load the image");
	}

	double *image = allocate(h,w);
	for(int i=0; i<h; i++) {
		for(int j=0; j<w; j++) {
			image[i*w+j] = img[i*w+j]/255.0;
		}
	}
	stbi_image_free(img);

	double *new = allocate(h,w);
	svd(h,w,image,k,new);

	unsigned char *out = malloc((size_t)h*w);
	for(int i=0; i<h; i++) {
		for(int j=0; j<w; j++) {
			double value = new[i*w+j];
			if(value<0) {
				value=0;
			}
			if(value>1) {
				value=1;
			}
			out[i*w+j] = (unsigned char)(value*255);
		}
	}

	stbi_write_png(nf,w,h,1,out,w);
	printf("New image saved\n");

	double err;
	for(int i=0; i<h*w; i++) {
		err += pow(image[i]-new[i],2);
	}
	printf("The Frobenius error is: %lf\n",sqrt(err));
	
	free(image);
	free(new);
	free(out);

	return 0;
}
