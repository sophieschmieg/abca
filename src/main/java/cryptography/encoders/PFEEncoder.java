package cryptography.encoders;

import java.math.BigInteger;

import cryptography.interfaces.Encoder;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers.IntE;
import util.MiscAlgorithms;

public class PFEEncoder implements Encoder<PFE> {
	private int byteLength;

	public PFEEncoder(PrimeField field) {
		this(field.characteristic());
	}

	public PFEEncoder(IntE prime) {
		this(prime.getValue());
	}

	public PFEEncoder(BigInteger prime) {
		this.byteLength = 8 * MiscAlgorithms.DivRoundUp(prime.bitLength(), 8);
	}

	public byte[] encode(PFE t) {
		byte[] result = new byte[byteLength];
		for (int i = 0; i < t.getValue().bitLength(); i++) {
			result[i / 8] |= t.getValue().testBit(i) ? 1 << (i % 8) : 0;
		}
		return result;
	}

}
