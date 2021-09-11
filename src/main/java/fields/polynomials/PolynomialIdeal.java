package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractIdeal;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;

public class PolynomialIdeal<T extends Element<T>> extends AbstractIdeal<Polynomial<T>> {
	private List<Polynomial<T>> basis;
	private PolynomialRing<T> polynomialRing;

	PolynomialIdeal(PolynomialRing<T> polynomialRing, List<Polynomial<T>> generators) {
		super(polynomialRing);
		this.polynomialRing = polynomialRing;
		this.basis = polynomialRing.buchberger(generators).getBasis();
	}

	@Override
	public boolean contains(Polynomial<T> polynomial) {
		return polynomialRing.reduce(polynomial, this.basis).equals(polynomialRing.zero());
	}

	@Override
	public Polynomial<T> residue(Polynomial<T> f) {
		return polynomialRing.reduce(f, this.basis);
	}

	public PolynomialIdeal<T> saturate(Polynomial<T> by) {
		PolynomialRing<T> ring = polynomialRing.addVariableWithElimination(1);
		Polynomial<T> t = ring.getVar(1);
		Polynomial<T> f = ring.getEmbeddingWithElimination(by, 1);
		Polynomial<T> inv = ring.subtract(ring.one(), ring.multiply(f, t));
		List<Polynomial<T>> generators = new ArrayList<Polynomial<T>>();
		generators.add(inv);
		for (Polynomial<T> p : this.basis) {
			generators.add(ring.getEmbeddingWithElimination(p, 1));
		}
		Ideal<Polynomial<T>> saturationIdeal = ring.getIdeal(generators);
		List<Polynomial<T>> saturationGenerators = new ArrayList<Polynomial<T>>();
		for (Polynomial<T> b : saturationIdeal.generators()) {
			if (b.leadingMonomial().exponents()[0] == 0)
				saturationGenerators.add(polynomialRing.getEmbeddingWithElimination(b, -1));
		}
		return new PolynomialIdeal<T>(polynomialRing, saturationGenerators);

	}

	public String toString() {
		return this.basis.toString();
	}

	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() {
		throw new InfinityException();
	}

	@Override
	public List<Polynomial<T>> generators() {
		return Collections.unmodifiableList(this.basis);
	}

	@Override
	public List<Polynomial<T>> generate(Polynomial<T> t) {
		return polynomialRing.generalQuotientAndRemainder(t, basis).getQuotients();
	}

	public Ideal<T> intersectToRing() {
		List<T> ringGenerators = new ArrayList<>();
		for (Polynomial<T> generator : generators()) {
			if (generator.degree() == 0) {
				ringGenerators.add(generator.leadingCoefficient());
			}
		}
		return polynomialRing.getRing().getIdeal(ringGenerators);
	}

	@Override
	public boolean isPrime() {
		for (Polynomial<T> polynomial : generators()) {
			if (!polynomialRing.isIrreducible(polynomial)) {
				return false;
			}
		}
		return intersectToRing().isPrime();
	}

	@Override
	public boolean isMaximal() {
		return isPrime() && generators().size() == polynomialRing.numberOfVariables() && intersectToRing().isMaximal();
	}
	
	public boolean isMaximalOverAlgebraicClosure() {
		if (polynomialRing.numberOfVariables() != polynomialRing.krullDimension()) {
			throw new ArithmeticException("Only defined over field!");
		}
		if (polynomialRing.numberOfVariables() != generators().size()) {
			return false;
		}
		for (Polynomial<T> generator : generators()) {
			if (generator.degree() != 1) {
				return false;
			}
		}
		return true;
	}

}