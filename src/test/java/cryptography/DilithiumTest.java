package cryptography;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.security.SecureRandom;

import org.junit.jupiter.api.Test;

import cryptography.Dilithium.DilithiumPrivateKey;
import cryptography.Dilithium.DilithiumPublicKey;
import cryptography.Dilithium.DilithiumSignature;

class DilithiumTest {

	@Test
	void test() {
		Dilithium d3 = Dilithium.DILITHIUM3;
		DilithiumPrivateKey privateKey = d3.createPrivateKey();
		DilithiumPublicKey publicKey = d3.createPublicKey(privateKey);
		byte[] message = new byte[100];
		new SecureRandom().nextBytes(message);
		DilithiumSignature signature = d3.sign(message, privateKey);
		assertTrue(d3.verify(message, signature, publicKey));
		assertFalse(d3.verify(new byte[100], signature, publicKey));
	}

}
