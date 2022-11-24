package cryptography;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import cryptography.RsaKem.RsaPrivateKey;
import cryptography.RsaKem.RsaPublicKey;
import fields.integers.Integers.IntE;
import util.Pair;

class RsaKemTest {

	@Test
	void test() {
		RsaKem rsa = new RsaKem(2048, Sha3.SHAKE_256);
		for (int keys = 0; keys < 10; keys++) {
			RsaPrivateKey privateKey = rsa.createPrivateKey();
			RsaPublicKey publicKey = rsa.createPublicKey(privateKey);
			for (int tc = 0; tc < 20; tc++) {
				Pair<VariableLengthKey, IntE> encapsulated = rsa.encapsulate(publicKey);
				VariableLengthKey decapsulated = rsa.decapsulate(encapsulated.getSecond(), privateKey);
				assertEquals(encapsulated.getFirst(), decapsulated);
			}
		}
	}

}
