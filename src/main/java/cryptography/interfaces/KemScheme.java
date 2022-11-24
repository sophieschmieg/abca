package cryptography.interfaces;

import cryptography.VariableLengthKey;
import fields.interfaces.Element;
import util.Pair;

public interface KemScheme<K extends Element<K>, PK extends Element<PK>, SK extends Element<SK>, S extends KemScheme<K, PK, SK, S>>
		extends Scheme<PK, SK, S> {
	SK createPrivateKey();

	Pair<VariableLengthKey, K> encapsulate(PK publicKey);

	VariableLengthKey decapsulate(K kem, SK privateKey);
}
