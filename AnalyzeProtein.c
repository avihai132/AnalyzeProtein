#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MIN_ARG_NUM 2
#define BASE_ARG_IDX 1

#define MAX_ATOMS 20000
#define COORD_NUM 3
#define MAX_LINE_LENGTH 81      // 80 + 1 because fgets adds an extra null char
#define COORD_STR_SIZE 9        // 8 + 1 for null char
#define X_COORD_OFFSET 30       // 31 - 1 because char array starts from 0
#define Y_COORD_OFFSET 38       // 39 - 1
#define Z_COORD_OFFSET 46       // 47 - 1

#define X_COORD_IDX 0
#define Y_COORD_IDX 1
#define Z_COORD_IDX 2


#define INVALID_USAGE_MSG "Invlid usage. Please enter at least one file path." // todo: change msg
#define ATOM_KEYWORD "ATOM"

void analyzeFile(char *filePath);
int loadCoordinates(char *filePath, float coordinatesArr[][COORD_NUM]);
int extractCoordinates(float coordinates[], char line[]);
void calcCg(float coordinatesArr[][COORD_NUM], int atomsNum, float cg[]);
float calcRg(float coordinatesArr[][COORD_NUM], int atomsNum, float cg[]);
float calcMaxDistance(float coordinatesArr[][COORD_NUM], int atomsNum);
float getDistance(float point1[], float point2[]);


/**
 * The entry point of the program.
 * @param argc The program argument count.
 * @param argv The program argument values.
 * @return Exit code. todo: check the return value?
 */
int main(int argc, char *argv[])
{
    if (argc < MIN_ARG_NUM)
    {
        printf(INVALID_USAGE_MSG);
        return -1; // todo: exit?
    }

    for (int i = BASE_ARG_IDX; i < argc; i++)
    {
        analyzeFile(argv[i]);
    }

    return 0;
}

/**
 * Analyzes the given file. todo: edit this doc
 * @param filePath The path to the file to analyze.
 */
void analyzeFile(char *filePath)
{
    float coordinatesArr[MAX_ATOMS][COORD_NUM] = {{0}};
    float cg[COORD_NUM] = {0.0f};
    float rg = 0.0f;
    float dMax = 0.0f;
    int atomsNum = 0;

    atomsNum = loadCoordinates(filePath, coordinatesArr);
    calcCg(coordinatesArr, atomsNum, cg);
    rg = calcRg(coordinatesArr, atomsNum, cg);
    dMax = calcMaxDistance(coordinatesArr, atomsNum);

    printf("PDB file %s, %d atoms were read\n", filePath, atomsNum);
    printf("Cg = %.3f %.3f %.3f\n", cg[X_COORD_IDX], cg[Y_COORD_IDX], cg[Z_COORD_IDX]);
    printf("Rg = %.3f\n", rg);
    printf("Dmax = %.3f\n", dMax);
}

/**
 * Loads the coordinates from the given file to the given array.
 * @param filePath The file to load the coordinates from.
 * @param coordinatesArr The array to load the coordinates to.
 */
int loadCoordinates(char *filePath, float coordinatesArr[][COORD_NUM])
{
    FILE *fp;
    char line[MAX_LINE_LENGTH]; // buffer to read a line from the file
    int curCoordinate = 0; // the current coordinate that we are at in the coordinates array

    fp = fopen(filePath, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Unable to open the file."); // todo: define the msg as a const?
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL)
    {
        if (extractCoordinates(coordinatesArr[curCoordinate], line) == 0)
        {
            curCoordinate++;
        }
    }

    fclose(fp); // todo: maybe need to check return value and see if there is an error
    return curCoordinate;
}

/**
 * Extracts the x, y and z coordinates from the given line and adds them to the coordinates array.
 * @param coordinates The array to load the extracted coordinates to.
 * @param line The line to extract the coordinates from.
 */
int extractCoordinates(float coordinates[], char line[]) {
    char atomSubStr[5] = {0}; // todo: change 5 to constant?
    char coordBuffer[COORD_STR_SIZE] = {0};
    memcpy(atomSubStr, line, 4); // todo: change 4 to constant?
    atomSubStr[4] = '\0';

    if (strcmp(atomSubStr, ATOM_KEYWORD) == 0)
    {
        // valid ATOM line
        // extract X coord
        memcpy(coordBuffer, &line[X_COORD_OFFSET], COORD_STR_SIZE - 1);
        coordBuffer[COORD_STR_SIZE - 1] = '\0';
        coordinates[0] = strtof(coordBuffer, NULL);

        // extract Y coord
        memcpy(coordBuffer, &line[Y_COORD_OFFSET], COORD_STR_SIZE - 1);
        coordBuffer[COORD_STR_SIZE - 1] = '\0';
        coordinates[1] = strtof(coordBuffer, NULL);

        // extract Z coord
        memcpy(coordBuffer, &line[Z_COORD_OFFSET], COORD_STR_SIZE - 1);
        coordBuffer[COORD_STR_SIZE - 1] = '\0';
        coordinates[2] = strtof(coordBuffer, NULL);

        return 0;
    }

    return -1;
}

/**
 * Calculates the center of gravity of the protein.
 * The center of gravity is defined to be the average by each coordinate of all the atoms.
 * @param coordinatesArr The coordinates of the atoms.
 * @param atomsNum The number of atoms.
 * @param cg The variable that will hold the result.
 */
void calcCg(float coordinatesArr[][COORD_NUM], int atomsNum, float cg[])
{
    float sumX = 0.0f;
    float sumY = 0.0f;
    float sumZ = 0.0f;
    for (int i = 0; i < atomsNum; i++)
    {
        sumX += coordinatesArr[i][X_COORD_IDX];
        sumY += coordinatesArr[i][Y_COORD_IDX];
        sumZ += coordinatesArr[i][Z_COORD_IDX];
    }

    cg[X_COORD_IDX] = sumX / atomsNum;
    cg[Y_COORD_IDX] = sumY / atomsNum;
    cg[Z_COORD_IDX] = sumZ / atomsNum;
}

/**
 * Calculates the radius of gyration.
 * The radius of gyration is defined to be the square root of the average of the quadratic distance
 * from the center of gravity.
 * @param coordinatesArr The coordinates of the atoms.
 * @param atomsNum The number of atoms.
 * @param cg The center of gravity.
 * @return The radius of gyration.
 */
float calcRg(float coordinatesArr[][COORD_NUM], int atomsNum, float cg[])
{
    float sum = 0.0f;

    for (int i = 0; i < atomsNum; i++)
    {
        sum += powf(getDistance(cg, coordinatesArr[i]), 2.0f);
    }

    return sqrtf(sum / atomsNum);
}

float calcMaxDistance(float coordinatesArr[][COORD_NUM], int atomsNum)
{
    float maxDistance = 0.0f;
    float curDistance = 0.0f;

    for (int i = 0; i < atomsNum; i++)
    {
        for (int j = i + 1; j < atomsNum; j++)
        {
            curDistance = getDistance(coordinatesArr[i], coordinatesArr[j]);
            if (maxDistance < curDistance)
            {
                maxDistance = curDistance;
            }
        }
    }

    return maxDistance;
}

/**
 * Calculates the distance between two points with (x,y,z) coordinates.
 * @param point1 The first point.
 * @param point2 The second point.
 * @return The distance between the two points.
 */
float getDistance(float point1[], float point2[])
{
    float xDist = powf(point1[X_COORD_IDX] - point2[X_COORD_IDX], 2.0f);
    float yDist = powf(point1[Y_COORD_IDX] - point2[Y_COORD_IDX], 2.0f);
    float zDist = powf(point1[Z_COORD_IDX] - point2[Z_COORD_IDX], 2.0f);
    return sqrtf(xDist + yDist + zDist);
}