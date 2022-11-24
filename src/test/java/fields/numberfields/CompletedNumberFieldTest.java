package fields.numberfields;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.Value;
import fields.numberfields.CompletedNumberField.Ext;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.Vector;
import util.Pair;

class CompletedNumberFieldTest {

	@Test
	void testCyclomatic8() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField
				.getNumberField(rationalPolynomials.add(rationalPolynomials.getVarPower(8), rationalPolynomials.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		int[] primes = new int[] { 2, 3, 5, 7, 41, 43 };
		int[] accuracy = new int[] { 32, 8, 8, 4, 4, 4 };
		int[] round = new int[] { 32, 4, 4, 4, 2, 2 };
		List<NFE> integers = new ArrayList<>();
		integers.add(order.zero());
		integers.add(order.one());
		integers.add(order.getInteger(-1));
		integers.add(order.getInteger(2));
		integers.add(order.getInteger(-2));
		integers.addAll(order.getModuleGenerators());
		for (int i = 0; i < primes.length; i++) {
			int prime = primes[i];
			System.out.println(prime);
			for (NumberFieldIdeal primeIdeal : order.idealsOver(prime)) {
				System.out.println(primeIdeal);
				CompletedNumberField complete = order.localizeAndQuotient(primeIdeal).complete(accuracy[i]).getField();
				for (NFE integer : integers) {
					System.out.println(integer);
					Ext ext = complete.fromElement(integer);
					System.out.println(ext);
					assertEquals(integer, complete.roundToInteger(ext, round[i]));
				}
			}
		}
	}

	@Test
	void testSqrtMinus5() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(
				rationalPolynomials.add(rationalPolynomials.getVarPower(2), rationalPolynomials.getInteger(5)));
		NumberFieldIntegers order = nf.maximalOrder();
		int[] primes = new int[] { 2, 3, 5, 7, 11 };
		List<NFE> integers = new ArrayList<>();
		for (int i = -3; i < 4; i++) {
			for (int j = -3; j < 4; j++) {
				integers.add(order.fromVector(new Vector<>(z.getInteger(i), z.getInteger(j))));
			}
		}
		for (int prime : primes) {
			for (NumberFieldIdeal primeIdeal : order.idealsOver(prime)) {
				CompletedNumberField complete = order.localizeAndQuotient(primeIdeal).complete(16).getField();
				System.out.println(primeIdeal);
				for (NFE integer : integers) {
					Ext completed = complete.fromElement(integer);
					System.out.println(primeIdeal);
					System.out.println("Target: " + integer);
					System.out.println("Required accuracy: " + roundToInteger(complete, completed, integer));
					NFE maxDenominator = null;
					int maxAccuracy = -1;
					for (NFE denominator : integers) {
						Ext completedDenominator = complete.fromElement(denominator);
						if (complete.valuation(completedDenominator).compareTo(complete.valuation(completed)) > 0) {
							continue;
						}
						Ext fraction = complete.divide(completed, completedDenominator);
						int accuracy = roundToRational(complete, fraction, order, integer, denominator);
						if (accuracy > maxAccuracy) {
							maxAccuracy = accuracy;
							maxDenominator = denominator;
						}
					}
					System.out.println("Target:  " + integer + "/" + maxDenominator);
					System.out.println("Required accuracy: " + maxAccuracy);
				}
			}
		}
	}

	private int roundToInteger(CompletedNumberField complete, Ext completed, NFE target) {
		for (int accuracy = 1; true; accuracy++) {
			NFE rounded = complete.roundToInteger(completed, accuracy);
			assertTrue(complete.valuation(complete.subtract(completed, complete.fromElement(rounded)))
					.compareTo(new Value(accuracy)) >= 0);
			if (rounded.equals(target)) {
				return accuracy;
			}
		}
	}

	private int roundToRational(CompletedNumberField complete, Ext completed, NumberFieldIntegers order,
			NFE targetNumerator, NFE targetDenominator) {
		for (int accuracy = 1; true; accuracy++) {
			Pair<NFE, NFE> rounded = complete.roundToRational(completed, accuracy);
			Value v = complete.valuation(
					complete.subtract(complete.multiply(completed, complete.fromElement(rounded.getSecond())),
							complete.fromElement(rounded.getFirst())));
			assertTrue(v.compareTo(new Value(accuracy)) >= 0);
			if (order.multiply(rounded.getFirst(), targetDenominator)
					.equals(order.multiply(rounded.getSecond(), targetNumerator))) {
				return accuracy;
			}
		}
	}
}
