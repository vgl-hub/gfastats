#include "zlib/zlib.h"
#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[])
{
	const unsigned char pData[] = {
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220,
		123, 34, 37, 83, 24, 2, 98, 178, 57, 220
	};
	unsigned long nDataSize = 100;

	printf("Initial size: %d\n", nDataSize);

	unsigned long nCompressedDataSize = nDataSize;
	unsigned char * pCompressedData = new unsigned char[nCompressedDataSize];
	
	int nResult = compress2(pCompressedData, &nCompressedDataSize, pData, nDataSize, 9);

	if (nResult == Z_OK)
	{
		printf("Compressed size: %d\n", nCompressedDataSize);

		unsigned char * pUncompressedData = new unsigned char[nDataSize];
		nResult = uncompress(pUncompressedData, &nDataSize, pCompressedData, nCompressedDataSize);
		if (nResult == Z_OK)
		{
			printf("Uncompressed size: %d\n", nDataSize);
			if (memcmp(pUncompressedData, pData, nDataSize) == 0)
				printf("Great Success\n");
		}
		delete [] pUncompressedData;
	}

	delete [] pCompressedData;

	return 0;
}