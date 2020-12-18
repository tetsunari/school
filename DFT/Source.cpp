#include<stdio.h>
#include<math.h>

double t_f[200001];
double zikai_f[200001];
//�t�[���G�ϊ�//
double F[101];
//��ꍀA[k],���B[k]//
double A[101];
double B[101];

void main() {

	///�e��p�����[�^///
	int k, n;
	int N = 200001;
	double  pi = 4.0 * atan(1.0);
	double dt = 5.0e-12;
	double T = N*dt;
	double df = 1.0 / T;

	//�t�@�C���ǂݍ��݂̍ۂ̕ϐ�//
	
	FILE *fp_r, *fp_w;

	//�t�@�C���̓ǂݍ���//
	fopen_s(&fp_r, "���E�̎��ԕω�(2-3�� 0.4MHz).txt", "r");
	for (n = 0; n <= 200000; n++) {

		fscanf_s(fp_r, "%lf %lf", &t_f[n], &zikai_f[n]);

		printf("%lf\n", zikai_f[n]);
	}
	fclose(fp_r);

	//DFT//
	fopen_s(&fp_w, "�t�[���G�ϊ�(2-3�� 0.4MHz).txt", "w");
	for (k = 0; k < 100; k++) {
		A[k] = 0.0;
		B[k] = 0.0;
		for (n = 0; n <200000; n++) {
			A[k] += (2.0 / N)*zikai_f[n] * cos(2.0*pi*n*k / N);
			B[k] += (2.0 / N)*zikai_f[n] * sin(2.0*pi*n*k / N);
		}
		F[k] = 20.0*log10(sqrt(pow(A[k], 2.0) + pow(B[k], 2.0)));
		//F[k] = sqrt(pow(A[k], 2.0) + pow(B[k], 2.0));
		fprintf(fp_w, "%e\t %e\n", k*df,F[k]);
		printf("%e %e\n", k*df, F[k]);
	}
	fclose(fp_w);
}