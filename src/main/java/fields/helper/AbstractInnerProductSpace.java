package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.MathMap;
import fields.interfaces.Ring;
import fields.interfaces.ValueField;

public abstract class AbstractInnerProductSpace<T extends Element<T>, S extends Element<S>> extends AbstractModule<T, S>
		implements InnerProductSpace<T, S> {
	private static Reals r = Reals.r(1024);

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	public Ring<T> getRing() {
		return getValueField();
	}

	public Field<T> getField() {
		return getValueField();
	}

	public List<S> getGenerators() {
		return getBasis();
	}

	public static <T extends Element<T>, S extends Element<S>> List<S> latticeReduction(InnerProductSpace<T, S> space,
			MathMap<T, T> round, List<S> s) {
		ValueField<T> f = space.getValueField();
		Real delta = r.getDouble(0.75);
		Real half = r.getDouble(0.5);
		List<S> basis = new ArrayList<>();
		basis.addAll(s);
		List<S> orthogonal = space.gramSchmidt(basis);
		List<List<T>> coefficients = new ArrayList<>();
		for (int i = 0; i < basis.size(); i++) {
			coefficients.add(new ArrayList<>());
			for (int j = 0; j < basis.size(); j++) {
				coefficients.get(i).add(f.zero());
			}
		}
		computeLLLCoefficients(space, basis, orthogonal, coefficients);
		int k = 1;
		while (k < basis.size()) {
			for (int j = k - 1; j >= 0; j--) {
				if (f.value(coefficients.get(k).get(j)).compareTo(half) > 0) {
					basis.set(k, space.subtract(basis.get(k),
							space.scalarMultiply(round.evaluate(coefficients.get(k).get(j)), basis.get(j))));
					orthogonal = space.gramSchmidt(basis);
					computeLLLCoefficients(space, basis, orthogonal, coefficients);
				}
			}
			Real c = f.value(coefficients.get(k).get(k - 1));
			Real ksquare = f.value(space.innerProduct(orthogonal.get(k), orthogonal.get(k)));
			Real km1square = f.value(space.innerProduct(orthogonal.get(k - 1), orthogonal.get(k - 1)));
			if (ksquare.compareTo(r.multiply(r.subtract(delta, r.multiply(c, c)), km1square)) >= 0) {
				k++;
			} else {
				S tmp = basis.get(k);
				basis.set(k, basis.get(k - 1));
				basis.set(k - 1, tmp);
				orthogonal = space.gramSchmidt(basis);
				computeLLLCoefficients(space, basis, orthogonal, coefficients);
				k = Math.max(k - 1, 1);
			}
		}
		return basis;
	}

	private static <T extends Element<T>, S extends Element<S>> void computeLLLCoefficients(
			InnerProductSpace<T, S> space, List<S> basis, List<S> orthogonal, List<List<T>> coefficients) {
		ValueField<T> f = space.getValueField();
		for (int i = 0; i < basis.size(); i++) {
			for (int j = 0; j < basis.size(); j++) {
				coefficients.get(i).set(j, f.divide(space.innerProduct(basis.get(i), orthogonal.get(j)),
						space.innerProduct(orthogonal.get(j), orthogonal.get(j))));
			}
		}
	}

}
