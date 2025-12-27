#include <mpi.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Function to print an array
void printArray(int* arr, int size, const string& name) {
    cout << name << ": ";
    for (int i = 0; i < size; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;
}

// Function to convert character to integer (A=0, B=1, ..., Z=25)
int charToInt(char c) {
    return toupper(c) - 'A';
}

// Function to convert integer to character (0=A, 1=B, ..., 25=Z)
char intToChar(int n) {
    return (n % 26 + 26) % 26 + 'A';
}

// Extended Euclidean Algorithm for modular inverse
int extendedGCD(int a, int b, int &x, int &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    
    int x1, y1;
    int gcd = extendedGCD(b, a % b, x1, y1);
    
    x = y1;
    y = x1 - (a / b) * y1;
    
    return gcd;
}

// Function to calculate modular multiplicative inverse
int modInverse(int a, int m) {
    int x, y;
    int gcd = extendedGCD(a, m, x, y);
    
    if (gcd != 1) {
        return -1; // Inverse doesn't exist
    }
    
    return (x % m + m) % m;
}

// Function to calculate the determinant of a 2x2 matrix
int calculateDeterminant(int matrix[2][2]) {
    return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
}

// Function to calculate the modular inverse of the determinant
int calculateInvDeterminant(int det) {
    // Make determinant positive
    det = (det % 26 + 26) % 26;
    
    if (det == 0) {
        cerr << "Error: Determinant is 0, matrix is not invertible!" << endl;
        return 0;
    }
    
    int inv = modInverse(det, 26);
    if (inv == -1) {
        cerr << "Error: Determinant " << det << " is not invertible modulo 26!" << endl;
        cerr << "Determinant must be coprime with 26 (gcd(det, 26) = 1)" << endl;
        return 0;
    }
    
    return inv;
}

// Function to calculate the inverse of a 2x2 matrix mod 26
void calculateInverse(int matrix[2][2], int inverse[2][2]) {
    // Calculate determinant
    int det = calculateDeterminant(matrix);
    int det_mod = (det % 26 + 26) % 26;
    cout << "Determinant: " << det << " (mod 26 = " << det_mod << ")" << endl;
    
    // Calculate modular inverse of determinant
    int detInv = calculateInvDeterminant(det);
    
    if (detInv == 0) {
        cerr << "Error: Matrix is not invertible modulo 26!" << endl;
        return;
    }
    
    cout << "Inverse of determinant: " << detInv << endl;
    
    // Calculate the adjugate matrix
    // For matrix [a b; c d], adjugate is [d -b; -c a]
    int adj[2][2];
    adj[0][0] = matrix[1][1];
    adj[0][1] = (-matrix[0][1] + 26) % 26;
    adj[1][0] = (-matrix[1][0] + 26) % 26;
    adj[1][1] = matrix[0][0];
    
    cout << "\nAdjugate matrix (before multiplication):" << endl;
    cout << "[ " << adj[0][0] << " " << adj[0][1] << " ]" << endl;
    cout << "[ " << adj[1][0] << " " << adj[1][1] << " ]" << endl;
    
    // Multiply adjugate by inverse of determinant mod 26
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            inverse[i][j] = (adj[i][j] * detInv) % 26;
            // Ensure positive
            inverse[i][j] = (inverse[i][j] + 26) % 26;
        }
    }
}

// Function to decode a 2-integer block using matrix multiplication mod 26
void decodeBlock(int block[2], int inverseMatrix[2][2], int result[2]) {
    // Perform matrix multiplication: result = inverseMatrix * block mod 26
    result[0] = (inverseMatrix[0][0] * block[0] + inverseMatrix[0][1] * block[1]) % 26;
    result[1] = (inverseMatrix[1][0] * block[0] + inverseMatrix[1][1] * block[1]) % 26;
    
    // Ensure positive results
    result[0] = (result[0] + 26) % 26;
    result[1] = (result[1] + 26) % 26;
}

// Verify inverse matrix
void verifyInverse(int matrix[2][2], int inverse[2][2]) {
    cout << "\nVerifying K * K⁻¹ = I (mod 26):" << endl;
    
    int product[2][2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            product[i][j] = 0;
            for (int k = 0; k < 2; k++) {
                product[i][j] += matrix[i][k] * inverse[k][j];
            }
            product[i][j] = (product[i][j] % 26 + 26) % 26;
        }
    }
    
    cout << "[ " << product[0][0] << " " << product[0][1] << " ]" << endl;
    cout << "[ " << product[1][0] << " " << product[1][1] << " ]" << endl;
    
    if (product[0][0] == 1 && product[0][1] == 0 && 
        product[1][0] == 0 && product[1][1] == 1) {
        cout << "✓ Inverse matrix is correct!" << endl;
    } else {
        cout << "✗ Inverse matrix verification failed!" << endl;
        cout << "Expected identity matrix!" << endl;
    }
}

// Manual test of decryption
void testDecryptionManually() {
    cout << "\n=== MANUAL DECRYPTION TEST ===" << endl;
    
    // Test with known values
    int inverse[2][2] = {{15, 17}, {20, 9}};
    
    // First block: HX -> 7, 23
    int block1[2] = {7, 23};
    int result1[2];
    decodeBlock(block1, inverse, result1);
    
    cout << "Block HX (7,23) -> (" << result1[0] << "," << result1[1] << ") -> " 
         << intToChar(result1[0]) << intToChar(result1[1]) << endl;
    
    // Second block: KS -> 10, 18
    int block2[2] = {10, 18};
    int result2[2];
    decodeBlock(block2, inverse, result2);
    
    cout << "Block KS (10,18) -> (" << result2[0] << "," << result2[1] << ") -> " 
         << intToChar(result2[0]) << intToChar(result2[1]) << endl;
    
    // Expected: E(4) A(0) T(19) T(19) ...
    cout << "Expected first 4 chars: E(4) A(0) T(19) T(19)" << endl;
}

// Coordinator (Master) function
void coordinator(int world_size) {
    auto start_time = high_resolution_clock::now();
    
    cout << "========================================" << endl;
    cout << "   PARALLEL HILL CIPHER DECRYPTION" << endl;
    cout << "========================================\n" << endl;
    
    // Given ciphertext
    string ciphertext = "HXKSSOEQYUQRWJRO";
    int cipherLength = ciphertext.length();
    
    cout << "Ciphertext to decrypt: " << ciphertext << endl;
    cout << "Ciphertext length: " << cipherLength << " characters\n" << endl;
    
    // Get encryption matrix from user
    int encryptMatrix[2][2];
    cout << "Enter 2x2 Hill Cipher matrix (row-wise, space separated): ";
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            cin >> encryptMatrix[i][j];
        }
    }
    
    cout << "\nKey Matrix:" << endl;
    cout << "[ " << encryptMatrix[0][0] << " " << encryptMatrix[0][1] << " ]" << endl;
    cout << "[ " << encryptMatrix[1][0] << " " << encryptMatrix[1][1] << " ]" << endl;
    
    // Calculate inverse matrix
    int inverseMatrix[2][2];
    calculateInverse(encryptMatrix, inverseMatrix);
    
    cout << "\nInverse Matrix (mod 26):" << endl;
    cout << "[ " << inverseMatrix[0][0] << " " << inverseMatrix[0][1] << " ]" << endl;
    cout << "[ " << inverseMatrix[1][0] << " " << inverseMatrix[1][1] << " ]" << endl;
    
    // Verify inverse
    verifyInverse(encryptMatrix, inverseMatrix);
    
    // Manual test
    testDecryptionManually();
    
    // Convert ciphertext to integer array
    int* cipherInts = new int[cipherLength];
    for (int i = 0; i < cipherLength; i++) {
        cipherInts[i] = charToInt(ciphertext[i]);
    }
    
    cout << "\nCiphertext as integers: ";
    printArray(cipherInts, cipherLength, "");
    
    // Calculate number of blocks (each block has 2 integers)
    int numBlocks = cipherLength / 2;
    cout << "\nTotal blocks to process: " << numBlocks << endl;
    
    // Broadcast inverse matrix to all participants
    MPI_Bcast(inverseMatrix, 4, MPI_INT, 0, MPI_COMM_WORLD);
    cout << "Broadcasted inverse matrix to all workers.\n" << endl;
    
    // Distribute blocks to participants
    int numWorkers = world_size - 1;
    
    // Adjust for single worker
    int baseBlocks, remainder;
    if (numWorkers == 0) {
        numWorkers = 1; // At least pretend we have workers
        baseBlocks = numBlocks;
        remainder = 0;
    } else {
        baseBlocks = numBlocks / numWorkers;
        remainder = numBlocks % numWorkers;
    }
    
    cout << "Work Distribution:" << endl;
    cout << "Number of workers: " << numWorkers << endl;
    cout << "Base blocks per worker: " << baseBlocks << endl;
    
    // Send blocks to each participant
    int currentBlock = 0;
    for (int rank = 1; rank <= numWorkers; rank++) {
        // Calculate blocks for this participant
        int blocksForThisRank = baseBlocks;
        if (rank <= remainder) {
            blocksForThisRank++;
        }
        
        // Send number of blocks first
        if (numWorkers > 1) {
            MPI_Send(&blocksForThisRank, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
        }
        
        // Send the actual blocks
        if (blocksForThisRank > 0) {
            cout << "Sent " << blocksForThisRank << " blocks to worker " << rank 
                 << " (blocks " << currentBlock << "-" << currentBlock + blocksForThisRank - 1 << ")" << endl;
            
            if (numWorkers > 1) {
                for (int b = 0; b < blocksForThisRank; b++) {
                    int block[2] = {cipherInts[currentBlock*2], cipherInts[currentBlock*2 + 1]};
                    MPI_Send(block, 2, MPI_INT, rank, 1, MPI_COMM_WORLD);
                    currentBlock++;
                }
            } else {
                // Single process - process directly
                currentBlock += blocksForThisRank;
            }
        }
    }
    
    // Allocate array for decrypted text
    int* decryptedInts = new int[cipherLength];
    
    // For single process, decode directly
    if (numWorkers == 1) {
        for (int b = 0; b < numBlocks; b++) {
            int block[2] = {cipherInts[b*2], cipherInts[b*2 + 1]};
            int decodedBlock[2];
            decodeBlock(block, inverseMatrix, decodedBlock);
            decryptedInts[b*2] = decodedBlock[0];
            decryptedInts[b*2 + 1] = decodedBlock[1];
        }
        cout << "Worker processed " << numBlocks << " blocks (single process mode)." << endl;
    } else {
        // Receive decoded blocks from participants
        currentBlock = 0;
        for (int rank = 1; rank <= numWorkers; rank++) {
            int blocksForThisRank = baseBlocks;
            if (rank <= remainder) {
                blocksForThisRank++;
            }
            
            // Receive decoded blocks from this participant
            for (int b = 0; b < blocksForThisRank; b++) {
                int decodedBlock[2];
                MPI_Recv(decodedBlock, 2, MPI_INT, rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                // Store in result array
                decryptedInts[currentBlock*2] = decodedBlock[0];
                decryptedInts[currentBlock*2 + 1] = decodedBlock[1];
                currentBlock++;
            }
            
            if (blocksForThisRank > 0) {
                cout << "Worker " << rank << " processed " << blocksForThisRank << " blocks." << endl;
            }
        }
    }
    
    // Convert decrypted integers back to characters
    string plaintext = "";
    for (int i = 0; i < cipherLength; i++) {
        plaintext += intToChar(decryptedInts[i]);
    }
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
    
    // Output results
    cout << "\n==================================================" << endl;
    cout << "DECRYPTION RESULTS" << endl;
    cout << "==================================================" << endl;
    cout << "Ciphertext:  " << ciphertext << endl;
    cout << "Plaintext:   " << plaintext << endl;
    cout << "Decoded values: ";
    for (int i = 0; i < cipherLength; i++) {
        cout << decryptedInts[i] << " ";
    }
    cout << endl << endl;
    
    // Try to make it more readable
    string readable = plaintext;
    // Remove trailing 'X' characters (common padding in Hill Cipher)
    while (!readable.empty() && readable.back() == 'X') {
        readable.pop_back();
    }
    
    cout << "Readable form: " << readable << endl;
    cout << "==================================================\n" << endl;
    
    cout << "Execution time: " << duration.count() / 1000.0 << " ms" << endl;
    
    // Clean up
    delete[] cipherInts;
    delete[] decryptedInts;
}

// Participant (Worker) function
void participant(int world_rank) {
    // Receive inverse matrix from coordinator
    int inverseMatrix[2][2];
    MPI_Bcast(inverseMatrix, 4, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Receive number of blocks to decode
    int numBlocks;
    MPI_Recv(&numBlocks, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // Process each block
    for (int b = 0; b < numBlocks; b++) {
        int block[2];
        MPI_Recv(block, 2, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Decode the block
        int decodedBlock[2];
        decodeBlock(block, inverseMatrix, decodedBlock);
        
        // Send decoded block back to coordinator
        MPI_Send(decodedBlock, 2, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // Call appropriate function based on rank
    if (world_rank == 0) {
        coordinator(world_size);
    } else {
        participant(world_rank);
    }
    
    // Finalize MPI
    MPI_Finalize();
    
    return 0;
}