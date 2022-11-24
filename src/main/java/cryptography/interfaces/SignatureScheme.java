package cryptography.interfaces;

import fields.interfaces.Element;

public interface SignatureScheme<T extends Element<T>, PK extends Element<PK>, SK extends Element<SK>, S extends SignatureScheme<T, PK, SK, S>>
		extends Scheme<PK, SK, S> {
	SK createPrivateKey();
	T sign(byte[] message, SK privateKey);
	boolean verify(byte[] message, T signature, PK publicKey);
}
