package cryptography;

import org.junit.jupiter.api.Test;

import cryptography.interfaces.HashFunction;

class Sha3Test {

	@Test
	void test() {
		HashFunction sha3 = Sha3.SHA3_256;
		sha3.evaluate(new byte[] {-1, 15, 23, 25});
	}

}
