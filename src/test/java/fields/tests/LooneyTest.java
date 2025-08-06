package fields.tests;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Ring.ChineseRemainderPreparation;

class LooneyTest {

	private IntE randPrime(IntE min, IntE max) {
		Integers z = Integers.z();
		IntE p;
		do {
			p = z.add(min, z.getRandomElement(z.subtract(max, min)));
		} while (!z.isPrime(p));
		return p;
	}

	@Test
	void test() {
		Integers z = Integers.z();
		IntE max = z.getInteger(65536);
		IntE min = z.getInteger(64);
		IntE p = randPrime(min, max);
		IntE q = randPrime(min,  max);
		IntE n = z.multiply(p, q/*, r*/);
		List<IntE> moduli = new ArrayList<>();
		moduli.add(p);
		moduli.add(q);
		ChineseRemainderPreparation<IntE> prep = z.prepareChineseRemainderTheoremModuli(moduli);
		List<IntE> sqrtList = new ArrayList<>();
		sqrtList.add(z.one());
		sqrtList.add(z.negative(z.one()));
		IntE sqrt1 = z.chineseRemainderTheorem(sqrtList, prep);
		System.out.println("p = " + p);
		System.out.println("q = " + q);
		System.out.println("n = " + n);
		System.out.println("sqrt1 = " + sqrt1);
		System.out.println("sqrt1 * sqrt1 = " + z.remainder(z.multiply(sqrt1, sqrt1), n));
		IntE a = z.multiply(p, q);
		IntE b = sqrt1;
		PrimeField fp = PrimeField.getPrimeField(p);
		PrimeField fq = PrimeField.getPrimeField(q);
		while (!b.equals(z.zero())) {
	//		System.out.println(b + "^2 mod n = " + z.remainder(z.multiply(b, b), n));
	//		System.out.println(b + " mod p = " + z.centeredLift(fp.getInteger(b), p.getValue()));
	//		System.out.println(b + " mod q = " + z.centeredLift(fq.getInteger(b), q.getValue()));
			IntE remainder = z.remainder(a, b);
			IntE quotient = z.divideChecked(z.subtract(a, remainder), b);
			System.out.println(a + " = " + quotient + "*" + b + " + " + remainder);
				a = b;
			b = remainder;
		}
	}

}
