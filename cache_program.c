#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include<complex.h>
#include<stdbool.h>

//#define complex std::complex<double>

#define complex _Complex

//#define debug

#define MAX_LINES 65536
#define MAX_ASSOCIATIVITY 16

uint32_t cacheSize = (256*1024), addressBits = 32, cacheline;
uint8_t dataWidthBits = 32, dataWidthBytes = 4;

//read
uint64_t rmc, rbc, rbhc, rbmc, rbrdc, rbrc, flushCount, totalBytes_write=0, totalBytes_read=0;

//write
uint64_t wmc, wbc, wtc, wbhc, wbmc, wbrdc, wbrc;
#define WB  0
#define WTA 1
#define WTNA 2

//cache structure
typedef struct cache
{
	bool valid;
	bool dirty;
	uint8_t LRU;
	uint32_t tag;
	uint32_t data;
} cache;

cache cache_struct[MAX_LINES][MAX_ASSOCIATIVITY];


int get_Line(uint8_t *address, uint8_t burstLength, uint8_t N)
{
	uint32_t totalLines, lineBits, blockSize, blockBits, tagBits, line;
	blockSize = burstLength * dataWidthBytes;
	blockBits = log2(blockSize);
	totalLines = (cacheSize/(blockSize * N));
	lineBits = log2(totalLines);
	tagBits = addressBits - lineBits - blockBits;
	line = (((1<<(lineBits+blockBits)) - 1)&(uint32_t)address)>>blockBits;
	return line;
}

int get_Tag(uint8_t *address, uint8_t burstLength, uint8_t N)
{
	int tag = 0;
	uint32_t totalLines, lineBits, blockSize, blockBits, tagBits, l_value, totalTags;
	blockSize = burstLength * dataWidthBytes;
	blockBits = log2(blockSize);
	totalLines = (cacheSize/(  blockSize * N));
	lineBits = log2(totalLines);
	tagBits = addressBits - lineBits - blockBits;
	totalTags = pow(2, tagBits);
	tag = (int)address>>(lineBits+blockBits);
	return tag;

}

int cacheHit(uint32_t line,uint32_t tag, uint8_t N)
{
	int j, tagHit_index;
	for(j=0; j<N; j++)
	{
		if((cache_struct[line][j].tag == tag) && (cache_struct[line][j].valid)){
			tagHit_index = j;
			break;
		}
		else
			tagHit_index = -9999; //garbage value
	}
	return tagHit_index;
}


void clearCounters()
{
	wmc=0, wbc=0, wtc=0, wbhc=0, wbmc=0, wbrdc=0, wbrc=0;
	rmc=0, rbc=0, rbhc=0, rbmc=0, rbrdc=0, rbrc=0, flushCount=0;
}

void clearCache()
{
	int i = 0, j;
	for(i = 0; i<MAX_LINES; i++)
	{
		for(j=0;j<MAX_ASSOCIATIVITY; j++)
		{
			cache_struct[i][j].valid = false;
			cache_struct[i][j].dirty = false;
			cache_struct[i][j].LRU   =  j;
			cache_struct[i][j].tag = 0;
		}
	}
}

void flushCache()
{
	int i = 0, j;
	for(i = 0; i<MAX_LINES; i++)
	{
		for(j=0;j<MAX_ASSOCIATIVITY; j++)
		{
			if(cache_struct[i][j].valid && cache_struct[i][j].dirty == true)
			{
				flushCount++;
			}
		}
	}
}

int freeBlock(uint32_t line, uint8_t N)
{
	int i = 0;
	for(i = 0; i<N ; i++)
	{
		if(cache_struct[line][i].valid == 0)
			return i;
	}
	return -9999;
}

void updateLRU(uint32_t line, uint8_t N, uint8_t currentLRU)
{
	int i = 0;
	for(i = 0; i<N; i++)
	{
		if(cache_struct[line][i].LRU < currentLRU)
		{
			cache_struct[line][i].LRU++;
		}
	}
}


int lruBlock(uint32_t line, uint8_t N)
{
	int block_index = 0, i, max = -1;
	for(i = 0; i<N; i++)
	{
		if((max< (N-1)) && (cache_struct[line][i].LRU > max))
		{
			max = cache_struct[line][i].LRU;
			block_index = i;
			break;
		}
	}
	return block_index;
}



void readBlock(uint32_t line, uint32_t tag, uint8_t N)
{
	int rdBlock = 0, currentLRU = -1;
	rbc++;
	rdBlock = cacheHit(line,tag, N);
	if(rdBlock != -9999)
	{
		rbhc++;
		currentLRU = cache_struct[line][rdBlock].LRU;
		updateLRU(line, N, currentLRU);
		cache_struct[line][rdBlock].LRU = 0; //updating LRU of the current block to 0
	}
	else
	{
		rbmc++;
		rdBlock = freeBlock(line, N);
		if(rdBlock == -9999)
		{
			rdBlock = lruBlock(line, N);
			if(cache_struct[line][rdBlock].dirty == true)
			{
				rbrdc++;
				cache_struct[line][rdBlock].dirty = false;
			}
			rbrc++;
		}
		if(rdBlock != -9999)
		{
			cache_struct[line][rdBlock].valid = true;
			cache_struct[line][rdBlock].dirty = false;
			cache_struct[line][rdBlock].tag = tag;
			//update LRU
			currentLRU = cache_struct[line][rdBlock].LRU;
			updateLRU(line, N, currentLRU);
			cache_struct[line][rdBlock].LRU = 0;
		}
	}
}


void read_Memory(uint8_t *address, uint8_t byteCount, uint8_t burstLength, uint8_t N)
{
	int i;
	uint32_t old_line = -1;
	int line=0, tag=0;

	rmc++;

	for(i = 0; i<byteCount; i++)
	{
		line = get_Line(address, burstLength, N);
		tag = get_Tag(address, burstLength, N);


		if(line != old_line)
		{
			readBlock(line, tag, N);
			old_line = line;
		}
		address++;
	}
}

//write part

void writeBlock(uint32_t line, uint32_t tag, uint8_t N, uint8_t writeStrategy)
{
	int wrblock = 0, currentLRU;
	wbc++;
	if((writeStrategy==(WTA))||(writeStrategy==(WTNA)))
		wtc++;
	wrblock = cacheHit(line,tag, N);
	if(wrblock != -9999)
	{
		wbhc++;
	}
	else
	{
		wbmc++;
		wrblock = freeBlock(line, N);
		if((wrblock == -9999 ) && !(writeStrategy == WTNA))
		{
			wrblock = lruBlock(line, N);
			if(cache_struct[line][wrblock].dirty == true)
			{
				wbrdc++;
				cache_struct[line][wrblock].dirty = false;
			}
			wbrc++;
		}
	}
	if(wrblock != -9999)
	{
		currentLRU = cache_struct[line][wrblock].LRU;
		updateLRU(line, N, currentLRU);
		cache_struct[line][wrblock].LRU = 0;

		if(writeStrategy == WB)
		{
			cache_struct[line][wrblock].dirty = true;
		}

		cache_struct[line][wrblock].valid = true;
		cache_struct[line][wrblock].tag = tag;
	}
}


void write_Memory(uint8_t* address, uint8_t byteCount, uint8_t burstLength, uint8_t N, uint8_t writeStrategy)
{
	wmc++;
	int i =0;
	int old_line1 = -1, line, tag;
	for(i = 0; i < byteCount; i++)
	{
		line = get_Line(address, burstLength, N);
		tag = get_Tag(address, burstLength, N);

		if(line != old_line1)
		{
			writeBlock(line, tag, N, writeStrategy);
			old_line1 = line;
		}
		address++;
	}
}


void Radix2FFT_new(complex data[], int nPoints, int nBits, uint8_t BL, uint8_t N, uint8_t writeStrategy)
{
	// cooley-tukey radix-2, twiddle factor
	// adapted from Fortran code from Burrus, 1983
#pragma warning (disable: 4270)
	int i=0, j=0, k=0, l=0;
	int nPoints1, nPoints2;
	double complex cTemp, cTemp2;

	double dTheta, dDeltaCos,dDeltaSin, dCos, dSin;
#ifdef debug
	printf("Radix entry\n");
#endif
	nPoints2 = nPoints;
	read_Memory(&nPoints, sizeof(nPoints), BL, N);
	write_Memory(&nPoints2, sizeof(nPoints2), BL, N, writeStrategy);


	write_Memory(&k, sizeof(k),BL, N, writeStrategy);// for k=0 in the for loop
	read_Memory(&nBits, sizeof(nBits), BL, N);

	for (k = 1; k <= nBits; k++)
	{

		read_Memory(&k, sizeof(k), BL, N);//k<nbits & k++

		nPoints1 = nPoints2;
		read_Memory(&nPoints2, sizeof(nPoints2), BL, N);
		write_Memory(&nPoints1, sizeof(nPoints1), BL, N, writeStrategy);

		nPoints2 /= 2;
		read_Memory(&nPoints2, sizeof(nPoints2), BL, N);
		write_Memory(&nPoints2, sizeof(nPoints2), BL, N, writeStrategy);


		// Compute differential angles
		double dTheta = 2 * 3.14159257 / nPoints1;
		read_Memory(&nPoints1, sizeof(nPoints1), BL, N);
		write_Memory(&dTheta, sizeof(dTheta), BL, N, writeStrategy);

		double dDeltaCos = cos(dTheta);
		read_Memory(&dTheta, sizeof(dTheta), BL, N);
		write_Memory(&dDeltaCos, sizeof(dDeltaCos), BL, N, writeStrategy);


		double dDeltaSin = sin(dTheta);
		read_Memory(&dTheta, sizeof(dTheta), BL, N);
		write_Memory(&dDeltaSin, sizeof(dDeltaSin), BL, N, writeStrategy);


		// Initialize angles
		double dCos = 1;
		write_Memory(&dCos, sizeof(dCos), BL, N, writeStrategy);

		double dSin = 0;
		write_Memory(&dSin, sizeof(dSin), BL, N, writeStrategy);

		// Perform in-place FFT

		write_Memory(&j, sizeof(j), BL, N, writeStrategy);//j =0 next for loop
		read_Memory(&nPoints2, sizeof(nPoints2), BL, N);
		for (j = 0; j < nPoints2; j++)
		{


			read_Memory(&j, sizeof(j), BL, N);//j<npoints2 & j++

			i = j;
			read_Memory(&j, sizeof(j), BL, N);
			write_Memory(&i, sizeof(i), BL, N, writeStrategy);

			read_Memory(&i, sizeof(i), BL, N);
			read_Memory(&nPoints, sizeof(nPoints), BL, N);
			while (i < nPoints)
			{
				read_Memory(&i, sizeof(i), BL, N);

				read_Memory(&nPoints2, sizeof(nPoints2), BL, N);
				l = i + nPoints2;
				write_Memory(&l, sizeof(l), BL, N, writeStrategy);

				cTemp = data[i] - data[l];
				read_Memory(&i, sizeof(i), BL, N);
				read_Memory(&l, sizeof(l), BL, N);
				read_Memory(&data[i], sizeof(data[i]), BL, N);
				read_Memory(&data[l], sizeof(data[l]), BL, N);
				write_Memory(&cTemp, sizeof(cTemp), BL, N, writeStrategy);


				cTemp2 = data[i] + data[l];
				read_Memory(&i, sizeof(i), BL, N);
				read_Memory(&l, sizeof(l), BL, N);
				read_Memory(&data[i], sizeof(data[i]), BL, N);
				read_Memory(&data[l], sizeof(data[l]), BL, N);
				write_Memory(&cTemp2, sizeof(cTemp2), BL, N, writeStrategy);



				data[i] = cTemp2;
				read_Memory(&cTemp2, sizeof(cTemp2), BL, N);
				read_Memory(&i, sizeof(i), BL, N);
				write_Memory(&data[i], sizeof(data[i]), BL, N, writeStrategy);


				cTemp2 = CMPLX(dCos * creal(cTemp) + dSin * cimag(cTemp),dCos * cimag(cTemp) - dSin * creal(cTemp));
				read_Memory(&cTemp, sizeof(creal(cTemp)), BL, N);
				read_Memory(&dCos, sizeof(dCos), BL, N);
				read_Memory(&cTemp, sizeof(cimag(cTemp)), BL, N);
				read_Memory(&dSin, sizeof(dSin), BL, N);
				write_Memory(&cTemp2, sizeof(cTemp2), BL, N, writeStrategy);

				data[l] = cTemp2;
				read_Memory(&cTemp2, sizeof(cTemp2), BL, N);
				read_Memory(&l, sizeof(l), BL, N);
				write_Memory(&data[l], sizeof(data[l]), BL, N, writeStrategy);

				i += nPoints1;
				read_Memory(&nPoints1, sizeof(nPoints1), BL, N);
				read_Memory(&i, sizeof(i), BL, N);
				write_Memory(&i, sizeof(i), BL, N, writeStrategy);

				write_Memory(&i, sizeof(i), BL, N, writeStrategy); //i in the while loop
				read_Memory(&nPoints, sizeof(nPoints), BL, N);
			}

			double dTemp = dCos;
			read_Memory(&dCos, sizeof(dCos), BL, N);
			write_Memory(&dTemp, sizeof(dTemp), BL, N, writeStrategy);

			dCos = dCos * dDeltaCos - dSin * dDeltaSin;
			read_Memory(&dCos, sizeof(dCos), BL, N);
			read_Memory(&dDeltaCos, sizeof(dDeltaCos), BL, N);
			read_Memory(&dSin, sizeof(dSin), BL, N);
			read_Memory(&dDeltaSin, sizeof(dDeltaSin), BL, N);
			write_Memory(&dCos, sizeof(dCos), BL, N, writeStrategy);

			dSin = dTemp * dDeltaSin + dSin * dDeltaCos;
			read_Memory(&dTemp, sizeof(dTemp), BL, N);
			read_Memory(&dDeltaSin, sizeof(dDeltaSin), BL, N);
			read_Memory(&dSin, sizeof(dSin), BL, N);
			read_Memory(&dDeltaCos, sizeof(dDeltaCos), BL, N);
			write_Memory(&dSin, sizeof(dSin), BL, N, writeStrategy);

			read_Memory(&nPoints2, sizeof(nPoints2), BL, N);//j<npoints2
			write_Memory(&j, sizeof(j), BL, N, writeStrategy);//j++
		}


		read_Memory(&nBits, sizeof(nBits), BL, N);//k<nbits
		write_Memory(&k, sizeof(k), BL, N, writeStrategy);// k++
	}

	// Convert Bit Reverse Order to Normal Ordering
	j = 0;
	write_Memory(&j, sizeof(j), BL, N, writeStrategy);

	nPoints1 = nPoints - 1;
	read_Memory(&nPoints, sizeof(nPoints), BL, N);
	write_Memory(&nPoints1, sizeof(nPoints1), BL, N, writeStrategy);


	write_Memory(&i, sizeof(i), BL, N, writeStrategy);//for i=0
	read_Memory(&nPoints1, sizeof(nPoints1), BL, N);

	for (i = 0; i < nPoints1; i++)
	{

		read_Memory(&i, sizeof(i), BL, N);//i<npoints1 and i++
		read_Memory(&j, sizeof(j), BL, N);

		if (i < j) {

			cTemp = data[j];
			read_Memory(&j, sizeof(j), BL, N);
			read_Memory(&data[j], sizeof(data[j]), BL, N);
			write_Memory(&cTemp, sizeof(cTemp), BL, N, writeStrategy);

			cTemp2 = data[i];
			read_Memory(&i, sizeof(i), BL, N);
			read_Memory(&data[i], sizeof(data[i]), BL, N);
			write_Memory(&cTemp2, sizeof(cTemp2), BL, N, writeStrategy);

			data[i] = cTemp;
			read_Memory(&i, sizeof(i), BL, N);
			read_Memory(&cTemp, sizeof(cTemp), BL, N);
			write_Memory(&data[i], sizeof(data[i]), BL, N, writeStrategy);

			data[j] = cTemp2;
			read_Memory(&j, sizeof(j), BL, N);
			read_Memory(&cTemp2, sizeof(cTemp2), BL, N);
			write_Memory(&data[j], sizeof(data[j]), BL, N, writeStrategy);
		}

		k = nPoints / 2;
		read_Memory(&nPoints, sizeof(nPoints), BL, N);
		write_Memory(&k, sizeof(k), BL, N, writeStrategy);

		read_Memory(&j, sizeof(j), BL, N);
		read_Memory(&k, sizeof(k), BL, N);
		while (k <= j)
		{
			read_Memory(&k, sizeof(k), BL, N);

			j -= k;

			read_Memory(&j, sizeof(j), BL, N);
			write_Memory(&j, sizeof(j), BL, N, writeStrategy);


			k /= 2;
			read_Memory(&k, sizeof(k), BL, N);
			write_Memory(&k, sizeof(k), BL, N, writeStrategy);
			read_Memory(&j, sizeof(j), BL, N);
		}

		j += k;
		read_Memory(&k, sizeof(k), BL, N);
		read_Memory(&j, sizeof(j), BL, N);
		write_Memory(&j, sizeof(j), BL, N, writeStrategy);

		write_Memory(&i, sizeof(i), BL, N, writeStrategy);//for i++
		read_Memory(&nPoints1, sizeof(nPoints1), BL, N);
	}
#pragma warning(default: 4270)
}



int main(int argc, char* argv[])
{
	int i;
	double complex data[32768];


	FILE *fp;
	uint8_t BL, N, WS;

	int points = 32768;
    #define cycles 1


	int bits = ceil(log(points)/log(2));

	fp = fopen("/Users/sowmyasrinivas/Documents/DS/6313_cacheController/t16.csv", "w");
	if(fp == NULL)
	{
		printf("Failure");
		system("pause");
		return 1;
	}

	fprintf(fp, "BL,N,WS,rmc, rbc, rbhc, rbmc, rbrdc, rbrc,wmc, wbc, wtc, wbhc, wbmc, wbrdc, wbrc, flushCount\n");



	for(N=1; N<=16; N=N*2)
	{

		for(BL = 1; BL<=8; BL=BL*2)
		{

			for(WS=0; WS<=2; WS++)
			{
				clearCache();
				for (i = 0; i < points; i++)
				{
					data[i] = CMPLX(cos(2.0*3.1416*(float)i/(float)(points)*cycles), 0.0);
				}
				clearCounters();
				Radix2FFT_new(data, points, bits, BL, N, WS);
				printf("BL=%d\t N=%d\t WS=%d\n", BL, N, WS);
				flushCache();
				fprintf(fp, "%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n",BL,N,WS,rmc, rbc, rbhc, rbmc, rbrdc, rbrc, wmc, wbc, wtc, wbhc, wbmc,wbrdc,wbrc,flushCount);
			}
		}
	}

	printf("execution complete for 60 permutations\n");


	   for (i = 0; i < points; i++)
	   if(cimag(data[i]) >= 0.0)
	       printf("x[%d] = %2.4lf + j%2.4lf\n", i, creal(data[i]), cimag(data[i]));
	   else
	       printf("x[%d] = %2.4lf - j%2.4lf\n", i, creal(data[i]), - cimag(data[i]));
	return 0;
}

