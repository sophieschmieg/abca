package cryptography;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;

import org.junit.jupiter.api.Test;

import cryptography.Falcon.FalconPrivateKey;
import cryptography.Falcon.FalconPublicKey;
import cryptography.Falcon.FalconSignature;

class FalconTest {

	@Test
	void test() throws IOException {
		Falcon falcon = new Falcon(9);
		FalconPrivateKey privateKey = falcon.createPrivateKey();
		FalconPublicKey publicKey = falcon.createPublicKey(privateKey);
		FalconSignature signature = falcon.sign("I love my wife".getBytes(), privateKey);
		assertTrue(falcon.verify("I love my wife".getBytes(), signature, publicKey));
	}

}
