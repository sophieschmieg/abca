package cryptography.interfaces;

import fields.interfaces.Element;

public interface Scheme<PK extends Element<PK>, SK extends Element<SK>, S extends Scheme<PK, SK, S>> {
	PK createPublicKey(SK privateKey);
}
