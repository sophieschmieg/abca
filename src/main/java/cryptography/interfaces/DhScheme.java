package cryptography.interfaces;

import cryptography.VariableLengthKey;
import fields.interfaces.Element;

public interface DhScheme<PK extends Element<PK>, SK extends Element<SK>, S extends DhScheme<PK, SK, S>>
		extends Scheme<PK, SK, S> {
	public enum Role {
		Alice, Bob;
	}

	SK createPrivateKey(Role role);

	VariableLengthKey agree(SK privateKey, PK publicKey);
}
