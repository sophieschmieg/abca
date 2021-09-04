package fields.numberfields;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Ideal;
import fields.interfaces.Ring.FactorizationResult;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import util.Pair;

public class FractionalIdeal extends AbstractElement<FractionalIdeal> {
	private NumberFieldIntegers order;
	private NumberFieldIdeal numerator;
	private NumberFieldIdeal denominator;

	public FractionalIdeal(NumberFieldIntegers order, NumberFieldIdeal numerator) {
		this(order, numerator, order.getUnitIdeal());
	}
	
	public FractionalIdeal(NumberFieldIntegers order, NumberFieldIdeal numerator, IntE denominator) {
		this(order, numerator, order.getIdeal(Collections.singletonList(order.getInteger(denominator))));
	}

	public FractionalIdeal(NumberFieldIntegers order, NumberFieldIdeal numerator, NumberFieldIdeal denominator) {
		this.order = order;
		if (numerator.equals(order.getZeroIdeal())) {
			this.numerator = order.getZeroIdeal();
			this.denominator = order.getUnitIdeal();
			return;
		}
		FactorizationResult<Ideal<NFE>> numeratorFactors = order.idealFactorization(numerator);
		FactorizationResult<Ideal<NFE>> denominatorFactors = order.idealFactorization(denominator);
		Set<Ideal<NFE>> primeIdeals = new TreeSet<>();
		primeIdeals.addAll(numeratorFactors.primeFactors());
		primeIdeals.addAll(denominatorFactors.primeFactors());
		Map<NumberFieldIdeal, Integer> numeratorReduced = new TreeMap<>(); 
		Map<NumberFieldIdeal, Integer> denominatorReduced = new TreeMap<>(); 
		for (Ideal<NFE> prime : primeIdeals) {
			int multiplicity = numeratorFactors.multiplicity(prime) - denominatorFactors.multiplicity(prime);
			if (multiplicity > 0) {
				numeratorReduced.put((NumberFieldIdeal)prime, multiplicity);
			}
			if (multiplicity < 0) {
				denominatorReduced.put((NumberFieldIdeal)prime, -multiplicity);
			}
		}
		this.numerator = order.fromFactorization(numeratorReduced);
		this.denominator = order.fromFactorization(denominatorReduced);
	}
	
	@Override
	public String toString() {
		return getNumerator() + "/" + getDenominator();
	}
	
	public NumberFieldIdeal getNumerator() {
		return numerator;
	}
	
	public NumberFieldIdeal getDenominator() {
		return denominator;
	}

	@Override
	public int compareTo(FractionalIdeal o) {
		int cmp = denominator.compareTo(o.denominator);
		if (cmp != 0) {
			return cmp;
		}
		return numerator.compareTo(o.numerator);
	}

	public Pair<NumberFieldIdeal, IntE> clearDenominator() {
		Integers z = Integers.z();
		FactorizationResult<Ideal<NFE>> denominatorFactors = order.idealFactorization(denominator);
		Map<IntE, Integer> primes = new TreeMap<>();
		for (Ideal<NFE> primeFactor : denominatorFactors.primeFactors()) {
			NumberFieldIdeal ideal = (NumberFieldIdeal) primeFactor;
			if (!primes.containsKey(ideal.intGenerator())) {
				primes.put(ideal.intGenerator(), 0);
			}
			primes.put(ideal.intGenerator(),
					Math.max(primes.get(ideal.intGenerator()), denominatorFactors.multiplicity(primeFactor)));
		}
		IntE num = z.one();
		for (IntE prime : primes.keySet()) {
			num = z.multiply(num, z.power(prime, primes.get(prime)));
		}
		FractionalIdeal cleared = new FractionalIdeal(order,
				order.getIdeal(Collections.singletonList(order.getEmbedding(num))), denominator);
		if (!cleared.isInteger()) {
			throw new ArithmeticException("Could not clear the denominator");
		}
		NumberFieldIdeal numerator = order.multiply(this.numerator, cleared.numerator);
		return new Pair<>(numerator, num);
	}

	public boolean isPrincipal() {
		return clearDenominator().getFirst().isPrincipal();
	}

	public boolean isInteger() {
		return denominator.equals(order.getUnitIdeal());
	}
}
