package cryptography.interfaces;

import fields.interfaces.Element;
import util.Pair;

public interface KemScheme<T extends Element<T>, K extends Element<K>, PK extends Element<PK>, SK extends Element<SK>, S extends KemScheme<T, K, PK, SK, S>> extends Scheme<PK, SK, S> {
	SK createPrivateKey();
	Pair<T, K> encapsulate(PK publicKey);
	T decapsulate(K kem, SK privateKey);
}
