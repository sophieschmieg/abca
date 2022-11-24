package cryptography;

import cryptography.interfaces.SignatureScheme;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.NTTRing;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import util.Pair;

public class Dilithium implements SignatureScheme<ByteArray, ByteArray, ByteArray, Dilithium> {
	private IntE q;
	private PrimeField field;
	private int degree;
	private NTTRing<PFE, PrimeField> cylcotomic;

	public Dilithium() {
		this((1 << 23) - (1 << 13) + 1, 256, 1753);
	}

	public Dilithium(int q, int degree, int r) {
		Integers z = Integers.z();
		this.q = z.getInteger(q);
		this.degree = degree;
		this.field = PrimeField.getPrimeField(q);
		this.cylcotomic = new NTTRing<PrimeField.PFE, PrimeField>(field, degree, field.getInteger(r));
	}

	private Pair<IntE, IntE> power2Round(IntE t, int exp) {
		Integers z = Integers.z();
		t = z.remainder(t, q);
		IntE power = z.power(z.getInteger(2), exp);
		IntE halfPower = z.divideChecked(power, z.getInteger(2));
		IntE r0 = z.remainder(t, power);
		if (r0.compareTo(halfPower) > 0) {
			r0 = z.subtract(r0, power);
		}
		IntE quotient = z.divideChecked(z.subtract(t, r0), power);
		return new Pair<>(quotient, r0);
	}

	private Pair<IntE, IntE> decompose(IntE t, IntE alpha) {
		Integers z = Integers.z();
		t = z.remainder(t, q);
		IntE remainder = z.remainder(t, alpha);
		if (z.multiply(2, remainder).compareTo(alpha) > 0) {
			remainder = z.subtract(remainder, alpha);
		}
		if (z.subtract(t, remainder).equals(z.subtract(q, z.one()))) {
			return new Pair<>(z.zero(), z.subtract(remainder, z.one()));
		}
		IntE quotient = z.divideChecked(z.subtract(t, remainder), alpha);
		return new Pair<>(quotient, remainder);
	}

	private IntE highBits(IntE t, IntE alpha) {
		return decompose(t, alpha).getFirst();
	}

	private IntE lowBits(IntE t, IntE alpha) {
		return decompose(t, alpha).getSecond();
	}
	
	private boolean makeHint(IntE s, IntE t, IntE alpha) {
		IntE high = highBits(t, alpha);
	}
	
	private IntE useHint() {
		
	}
}
