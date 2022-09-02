package cryptography.interfaces;

public interface HashFunction {
	byte[] evaluate(byte[] input);
	int outputLength();
	int blockSize();
}
