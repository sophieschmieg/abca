package cryptography.encoders;

import cryptography.ByteArray;
import cryptography.interfaces.Encoder;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.FieldExtension;

public class FieldExtensionEncoder<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, FE extends FieldExtension<T, S, FE>>
		implements Encoder<S> {
	private int degree;
	private Encoder<T> encoder;

	public FieldExtensionEncoder(Encoder<T> encoder, FE field) {
		this.degree = field.degree();
		this.encoder = encoder;
	}

	public byte[] encode(S t) {
		ByteArray result = new ByteArray(new byte[0]);
		for (int i = 0; i < degree; i++) {
			result = ByteArray.concatenate(result.array(), encoder.encode(t.asPolynomial().univariateCoefficient(i)));
		}
		return result.array();
	}

}
