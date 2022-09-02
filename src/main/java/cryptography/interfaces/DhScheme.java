package cryptography.interfaces;

import fields.interfaces.Element;

public interface DhScheme<T extends Element<T>, PK extends Element<PK>, SK extends Element<SK>, S extends DhScheme<T, PK, SK, S>>
		extends Scheme<PK, SK, S> {
	public enum Role {
		Alice, Bob;
	}

	SK createPrivateKey(Role role);

	T agree(SK privateKey, PK publicKey);
}
