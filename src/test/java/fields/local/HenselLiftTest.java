package fields.local;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField.PAdicNumber;

class HenselLiftTest {

	@Test
	void test() {
		PAdicField z2 = new PAdicField(BigInteger.TWO, 32);
		UnivariatePolynomialRing<PAdicNumber> r = z2.getUnivariatePolynomialRing();
		UnivariatePolynomial<PAdicNumber> f = r.getPolynomial(z2.getInteger(2), z2.one(), z2.one());
		System.out.println(z2.factorization(f));
	}

}
