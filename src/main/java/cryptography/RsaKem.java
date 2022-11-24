package cryptography;

import java.util.ArrayList;
import java.util.List;

import cryptography.RsaKem.RsaPrivateKey;
import cryptography.RsaKem.RsaPublicKey;
import cryptography.interfaces.ExtendedOutputFunction;
import cryptography.interfaces.KemScheme;
import fields.finitefields.ModuloIntegerRing;
import fields.finitefields.PrimeField;
import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Ring.ChineseRemainderPreparation;
import util.Pair;

public class RsaKem implements KemScheme< IntE, RsaPublicKey, RsaPrivateKey, RsaKem> {
	public class RsaPublicKey extends AbstractElement<RsaPublicKey> {
		private IntE modulus;
		private IntE exponent;
		private ModuloIntegerRing mod;

		public RsaPublicKey(IntE modulus, IntE exponent) {
			this.modulus = modulus;
			this.exponent = exponent;
			mod = new ModuloIntegerRing(modulus.getValue());
		}

		public IntE getModulus() {
			return modulus;
		}

		public IntE getExponent() {
			return exponent;
		}

		@Override
		public int compareTo(RsaPublicKey o) {
			int cmp = modulus.compareTo(o.modulus);
			if (cmp != 0) {
				return cmp;
			}
			return exponent.compareTo(o.exponent);
		}
	}

	public class RsaPrivateKey extends AbstractElement<RsaPrivateKey> {
		private RsaPublicKey publicKey;
		private IntE eulerToitent;
		private IntE exponent;
		private ChineseRemainderPreparation<IntE> crt;
		private PrimeField mod1;
		private PrimeField mod2;

		private RsaPrivateKey(RsaPublicKey publicKey, IntE prime1, IntE prime2) {
			Integers z = Integers.z();
			this.publicKey = publicKey;
			mod1 = PrimeField.getPrimeField(prime1);
			mod2 = PrimeField.getPrimeField(prime2);
			if (!z.multiply(prime1, prime2).equals(publicKey.modulus)) {
				throw new RuntimeException("private key wrong!");
			}
			if (!z.isPrime(prime1) || !z.isPrime(prime2)) {
				throw new RuntimeException("private key wrong!");
			}
			if (prime1.equals(prime2)) {
				throw new RuntimeException("private key wrong!");
			}
			this.eulerToitent = z.multiply(z.subtract(prime1, z.one()), z.subtract(prime2, z.one()));
			ModuloIntegerRing modToitent = new ModuloIntegerRing(eulerToitent.getValue());
			this.exponent = modToitent.lift(modToitent.inverse(modToitent.reduce(publicKey.exponent)));
			List<IntE> primes = new ArrayList<>();
			primes.add(prime1);
			primes.add(prime2);
			this.crt = z.prepareChineseRemainderTheoremModuli(primes);
		}

		public RsaPublicKey getPublicKey() {
			return publicKey;
		}

		@Override
		public int compareTo(RsaPrivateKey o) {
			return publicKey.compareTo(o.publicKey);
		}
	}

	private int bitSize;
	private IntE primeMax;
	private IntE primeMin;
	private IntE diff;
	private ExtendedOutputFunction hash;

	public RsaKem(int bitSize, ExtendedOutputFunction hash) {
		if (bitSize % 2 != 0) {
			throw new RuntimeException("Expected even bit size!");
		}
		this.bitSize = bitSize;
		this.hash = hash;
		Integers z = Integers.z();
		primeMax = z.power(z.getInteger(2), bitSize / 2);
		primeMin = z.multiply(z.power(z.getInteger(2), bitSize / 2 - 2), z.getInteger(3));
		diff = z.subtract(primeMax, primeMin);
	}

	@Override
	public String toString() {
		return "RSA-" + bitSize + "-KEM-" + hash;
	}

	@Override
	public RsaPublicKey createPublicKey(RsaPrivateKey privateKey) {
		return privateKey.getPublicKey();
	}

	@Override
	public RsaPrivateKey createPrivateKey() {
		Integers z = Integers.z();
		IntE exponent = z.getInteger(65537);
		IntE prime1;
		do {
			prime1 = z.add(z.getRandomElement(diff), primeMin);
		} while (z.isDivisible(z.subtract(prime1, z.one()), exponent) || !z.isPrime(prime1));
		IntE prime2;
		do {
			prime2 = z.add(z.getRandomElement(diff), primeMin);
		} while (z.isDivisible(z.subtract(prime2, z.one()), exponent) || !z.isPrime(prime2));
		return new RsaPrivateKey(new RsaPublicKey(z.multiply(prime1, prime2), z.getInteger(65537)), prime1, prime2);
	}

	@Override
	public Pair<VariableLengthKey, IntE> encapsulate(RsaPublicKey publicKey) {
		Integers z = Integers.z();
		IntE secret = z.getRandomElement(publicKey.modulus);
		VariableLengthKey key = new VariableLengthKey(secret.getValue().toByteArray(), hash);
		IntE kem = publicKey.mod.lift(publicKey.mod.power(publicKey.mod.reduce(secret), publicKey.exponent));
		return new Pair<>(key, kem);
	}

	@Override
	public VariableLengthKey decapsulate(IntE kem, RsaPrivateKey privateKey) {
		List<IntE> powers = new ArrayList<>();
		powers.add(
				privateKey.mod1.liftToInteger(privateKey.mod1.power(privateKey.mod1.reduce(kem), privateKey.exponent)));
		powers.add(
				privateKey.mod2.liftToInteger(privateKey.mod2.power(privateKey.mod2.reduce(kem), privateKey.exponent)));
		IntE secret = Integers.z().chineseRemainderTheorem(powers, privateKey.crt);
		return new VariableLengthKey(secret.getValue().toByteArray(), hash);
	}
}
