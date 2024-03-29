package fields.vectors;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Lattice;
import fields.interfaces.RealInnerProductSpace;
import fields.interfaces.ValueField;
import fields.vectors.DualVectorSpace.Dual;
import util.MiscAlgorithms;

public abstract class AbstractRealInnerProductSpace<T extends Element<T>, S extends Element<S>>
		extends AbstractInnerProductSpace<T, S> implements RealInnerProductSpace<T, S> {

	@Override
	public T conjugateScalar(T t) {
		return t;
	}

	@Override
	public S conjugateVector(S t) {
		return t;
	}

	@Override
	public ComplexNumber asComplexNumber(T t) {
		return Complex.c(getValueField().getReals().precision()).getEmbedding(asReal(t));
	}

	@Override
	public T fromComplexNumber(ComplexNumber t) {
		return fromReal(t.realPart());
	}

	@Override
	public LinearProgramResult<T, S> linearProgram(Polytope<T, S> polytope, Dual<T, S> maximize) {
		return polytope.maximize(maximize);
	}

	@Override
	public <R extends Element<R>> List<R> latticeReduction(Lattice<R, T, S> lattice) {
		return latticeReduction(lattice.getModuleGenerators(), lattice);
	}

	@Override
	public <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice) {
		return latticeReduction(sublatticeBasis, lattice, true);
	}

	@Override
	public <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice,
			boolean isBasis) {
		if (sublatticeBasis.isEmpty()) {
			return Collections.emptyList();
		}
		return latticeReduction(sublatticeBasis, lattice, 0.75, 10, isBasis);
	}

	// @Override
	private <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice,
			double deltaAsDouble, boolean isBasis) {
		if (!lattice.getVectorSpace().equals(this)) {
			throw new ArithmeticException(lattice + " not in " + this + ", but in " + lattice.getVectorSpace());
		}
		List<R> prevResult;
		List<R> currentResult = sublatticeBasis;
		int counter = 0;
		do {
			counter++;
			prevResult = currentResult;
			currentResult = computeLatticeReduction(currentResult, lattice, deltaAsDouble, isBasis);
		} while (counter < 1000 && (prevResult == null || !prevResult.equals(currentResult)));
		return currentResult;
	}

	// @Override
	private <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice,
			double deltaAsDouble, int blockSize, boolean isBasis) {
		List<R> basis = latticeReduction(sublatticeBasis, lattice, deltaAsDouble, isBasis);
		ValueField<T> f = getValueField();
		List<S> orthogonal = new ArrayList<>();
		List<T> values = new ArrayList<>();
		List<List<T>> coefficients = new ArrayList<>();
		for (int i = 0; i < basis.size(); i++) {
			coefficients.add(new ArrayList<>());
			orthogonal.add(lattice.embedding(basis.get(i)));
			S embedded = orthogonal.get(i);
			for (int l = 0; l < i; l++) {
				coefficients.get(i).add(f.divide(innerProduct(embedded, orthogonal.get(l)), values.get(l)));
				orthogonal.set(i,
						subtract(orthogonal.get(i), scalarMultiply(coefficients.get(i).get(l), orthogonal.get(l))));
			}
			values.add(innerProduct(orthogonal.get(i), orthogonal.get(i)));
		}
		int z = 0;
		int j = 0;
		T eps = fromReal(f.getReals().getPowerOfTwo(-f.getReals().precision() + 1));
		while (z < basis.size() - 1) {
			j++;
			int k = Math.min(j + blockSize - 1, basis.size());
			if (j == basis.size()) {
				j = 1;
				k = Math.min(blockSize, basis.size());
			}
			R min = enumerate(j - 1, k - 1, basis, coefficients, values, lattice);
			S basisVector = lattice.embedding(basis.get(j - 1));
			S minVector = lattice.embedding(min);
			int h = Math.min(k, basis.size());
			List<R> subBasis;
			if (f.add(innerProduct(minVector, minVector), eps).compareTo(innerProduct(basisVector, basisVector)) < 0) {
				subBasis = new ArrayList<>();
				subBasis.addAll(basis.subList(0, h));
				subBasis.add(j - 1, min);
				subBasis = latticeReduction(subBasis, lattice, deltaAsDouble, false);
				z = 0;
			} else {
				z++;
				subBasis = latticeReduction(basis.subList(0, h), lattice, 1.0, isBasis);
				if (h == basis.size()) {
					return subBasis;
				}
			}
			for (int i = 0; i < subBasis.size(); i++) {
				basis.set(i, subBasis.get(i));
			}
			values.clear();
			orthogonal.clear();
			coefficients.clear();
			for (int i = 0; i < basis.size(); i++) {
				coefficients.add(new ArrayList<>());
				orthogonal.add(lattice.embedding(basis.get(i)));
				S embedded = orthogonal.get(i);
				for (int l = 0; l < i; l++) {
					coefficients.get(i).add(f.divide(innerProduct(embedded, orthogonal.get(l)), values.get(l)));
					orthogonal.set(i,
							subtract(orthogonal.get(i), scalarMultiply(coefficients.get(i).get(l), orthogonal.get(l))));
				}
				values.add(innerProduct(orthogonal.get(i), orthogonal.get(i)));
			}
		}
		return latticeReduction(basis, lattice, 1.0, isBasis);
	}

	private <R extends Element<R>> List<R> computeLatticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice,
			double deltaAsDouble, boolean isBasis) {
		ValueField<T> f = getValueField();
		Reals r = f.getReals();
		Real delta = r.getDouble(deltaAsDouble);
		Real eps = r.getPowerOfTwo(-r.precision() + 10);
		Real half = r.add(r.getDouble(0.5), eps);
		List<R> basis = new ArrayList<>();
		basis.addAll(sublatticeBasis);
		List<S> orthogonal = new ArrayList<>();
		List<T> values = new ArrayList<>();
		for (int i = 0; i < basis.size(); i++) {
			values.add(f.zero());
			orthogonal.add(zero());
		}
		List<List<T>> coefficients = new ArrayList<>();
		for (int i = 0; i < basis.size(); i++) {
			coefficients.add(new ArrayList<>());
			for (int j = 0; j < i; j++) {
				coefficients.get(i).add(f.zero());
			}
		}
		int k = 1;
		int kmax = 0;
		orthogonal.set(0, lattice.embedding(basis.get(0)));
		values.set(0, innerProduct(orthogonal.get(0), orthogonal.get(0)));
		while (k < basis.size()) {
			if ((!isBasis && asReal(values.get(k - 1)).compareTo(eps) < 0) || values.get(k - 1).equals(r.zero())) {
				basis.remove(k - 1);
				coefficients.remove(coefficients.size() - 1);
				orthogonal.remove(k - 1);
				values.remove(k - 1);
				k = Math.max(1, k - 1);
				kmax = Math.min(kmax, k - 1);
				continue;
			}
			if (k > kmax) {
				kmax = k;
				orthogonal.set(k, lattice.embedding(basis.get(k)));
				S embedded = orthogonal.get(k);
				for (int j = 0; j < k; j++) {
					coefficients.get(k).set(j, f.divide(innerProduct(embedded, orthogonal.get(j)), values.get(j)));
					orthogonal.set(k,
							subtract(orthogonal.get(k), scalarMultiply(coefficients.get(k).get(j), orthogonal.get(j))));
				}
				values.set(k, innerProduct(orthogonal.get(k), orthogonal.get(k)));
			}
			reduceLLLBasis(k, k - 1, coefficients, basis, lattice, half);
			if (!isBasis && asReal(values.get(k)).compareTo(eps) < 0) {
				basis.remove(k);
				coefficients.remove(coefficients.size() - 1);
				orthogonal.remove(k);
				values.remove(k);
				k = Math.max(1, k - 1);
				kmax = Math.min(kmax, k - 1);
				continue;
			}
			Real kValue = asReal(values.get(k));
			Real adjustedValue = r.multiply(
					r.subtract(delta,
							asReal(f.multiply(coefficients.get(k).get(k - 1), coefficients.get(k).get(k - 1)))),
					asReal(values.get(k - 1)));
			int maxExponent = Math.max(kValue.exponent() - r.precision(), adjustedValue.exponent() - r.precision()) + 2;
			if (r.abs(r.subtract(kValue, adjustedValue)).compareTo(r.getPowerOfTwo(maxExponent)) > 0
					&& kValue.compareTo(adjustedValue) < 0) {
				swapLLLBasis(k, coefficients, basis, orthogonal, values, kmax);
				k = Math.max(1, k - 1);
			} else {
				for (int j = k - 2; j >= 0; j--) {
					reduceLLLBasis(k, j, coefficients, basis, lattice, half);
				}
				k++;
			}
		}
		return basis;
	}

//	private <R extends Element<R>> void printDebug(List<R> basis, Lattice<R, T, S> lattice) {
//		System.out.println();
//		ValueField<T> f = getValueField();
//		Reals r =f.getReals();
//		Real half = r.divide(r.getInteger(1), r.getInteger(2));
//		List<S> orthogonal = gramSchmidt(embedList(basis, lattice));
//		for (int i = 0; i < orthogonal.size(); i++) {
//			System.out.println("Orthogonal(" + i + "): " + orthogonal.get(i));
//			System.out.println("Value(" + i + "): " + innerProduct(orthogonal.get(i), orthogonal.get(i)));
//			for (int j = 0; j < i; j++) {
//				T coeff = f.divide(innerProduct(lattice.embedding(basis.get(i)), orthogonal.get(j)),
//						innerProduct(orthogonal.get(j), orthogonal.get(j)));
//				System.out.println("Coefficient(" + i + "," + j + "): "
//						+ coeff);
//			if (r.abs(asReal(coeff)).compareTo(half) > 0) {
//				System.err.println(coeff + " is too large!");
//			}
//			}
//		}
//	}

	private <R extends Element<R>> void reduceLLLBasis(int k, int j, List<List<T>> coefficients, List<R> basis,
			Lattice<R, T, S> lattice, Real half) {
		ValueField<T> f = getValueField();
		Reals r = f.getReals();
		if (r.abs(asReal(coefficients.get(k).get(j))).compareTo(half) < 0) {
			return;
		}
		IntE m = round(coefficients.get(k).get(j));
		basis.set(k, lattice.subtract(basis.get(k), lattice.scalarMultiply(m, basis.get(j))));
		coefficients.get(k).set(j, f.subtract(coefficients.get(k).get(j), f.getInteger(m)));
		for (int i = 0; i < j; i++) {
			coefficients.get(k).set(i,
					f.subtract(coefficients.get(k).get(i), f.multiply(m, coefficients.get(j).get(i))));
		}
	}

	private <R extends Element<R>> void swapLLLBasis(int k, List<List<T>> coefficients, List<R> basis,
			List<S> orthogonal, List<T> values, int kmax) {
		ValueField<T> f = getValueField();
		R tmp = basis.get(k);
		basis.set(k, basis.get(k - 1));
		basis.set(k - 1, tmp);
		for (int j = 0; j < k - 1; j++) {
			T c = coefficients.get(k).get(j);
			coefficients.get(k).set(j, coefficients.get(k - 1).get(j));
			coefficients.get(k - 1).set(j, c);
		}
		T c = coefficients.get(k).get(k - 1);
		T value = f.add(values.get(k), f.multiply(c, c, values.get(k - 1)));
		coefficients.get(k).set(k - 1, f.divide(f.multiply(c, values.get(k - 1)), value));
		values.set(k, f.divide(f.multiply(values.get(k - 1), values.get(k)), value));
		values.set(k - 1, value);
		S tmpOrthogonal = orthogonal.get(k - 1);
		orthogonal.set(k - 1, add(orthogonal.get(k), scalarMultiply(c, orthogonal.get(k - 1))));
		orthogonal.set(k,
				subtract(tmpOrthogonal, scalarMultiply(coefficients.get(k).get(k - 1), orthogonal.get(k - 1))));
		for (int i = k + 1; i <= kmax; i++) {
			T m = coefficients.get(i).get(k);
			coefficients.get(i).set(k, f.subtract(coefficients.get(i).get(k - 1), f.multiply(m, c)));
			coefficients.get(i).set(k - 1,
					f.add(m, f.multiply(coefficients.get(k).get(k - 1), coefficients.get(i).get(k))));
		}
	}

	// https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.57.3331&rep=rep1&type=pdf
	private <R extends Element<R>> R enumerate(int j, int k, List<R> basis, List<List<T>> coefficients, List<T> values,
			Lattice<R, T, S> lattice) {
		Integers z = Integers.z();
		ValueField<T> field = getValueField();
		T minValue = values.get(j);
		List<T> minValues = new ArrayList<>(); // c
		minValues.add(field.zero());
		List<IntE> minParameters = new ArrayList<>(); // u
		minParameters.add(z.one());
		List<IntE> parameters = new ArrayList<>(); // tilde u
		parameters.add(z.one());
		List<T> minCoefficients = new ArrayList<>(); // y
		minCoefficients.add(field.zero());
		List<IntE> delta = new ArrayList<>(); // Delta
		delta.add(z.zero());
		List<IntE> sign = new ArrayList<>(); // delta
		sign.add(z.one());
		List<IntE> secondary = new ArrayList<>(); // v
		secondary.add(z.zero());
		int s = j;
		int t = j;
		for (int i = j + 1; i < k + 1; i++) {
			minValues.add(field.zero());
			minParameters.add(z.zero());
			parameters.add(z.zero());
			minCoefficients.add(field.zero());
			delta.add(z.zero());
			sign.add(z.one());
			secondary.add(z.zero());
		}
		minValues.add(field.zero());
		while (true) {
			T coeff = field.add(field.getInteger(parameters.get(t - j)), minCoefficients.get(t - j));
			minValues.set(t - j, field.add(minValues.get(t - j + 1), field.multiply(coeff, coeff, values.get(t - j))));
			if (minValues.get(t - j).compareTo(minValue) < 0) {
				if (t > j) {
					t--;
					T coefficient = field.zero();
					for (int i = t + 1; i <= s; i++) {
						coefficient = field.add(field.multiply(parameters.get(i - j), coefficients.get(i).get(t)),
								coefficient);
					}
					minCoefficients.set(t - j, coefficient);
					IntE rounded = asReal(field.negative(coefficient)).round();
					parameters.set(t - j, rounded);
					secondary.set(t - j, rounded);
					delta.set(t - j, z.zero());
					if (field.getInteger(rounded).compareTo(field.negative(coefficient)) > 0) {
						sign.set(t - j, z.getInteger(-1));
					} else {
						sign.set(t - j, z.one());
					}
				} else {
					minValue = minValues.get(0);
					minParameters.clear();
					for (IntE parameter : parameters) {
						minParameters.add(parameter);
					}
				}
			} else {
				t++;
				if (t > k) {
					break;
				}
				s = Math.max(s, t);
				if (t < s) {
					delta.set(t - j, z.negative(delta.get(t - j)));
				}
				if (z.multiply(sign.get(t - j), delta.get(t - j)).compareTo(z.zero()) >= 0) {
					delta.set(t - j, z.add(delta.get(t - j), sign.get(t - j)));
				}
				parameters.set(t - j, z.add(secondary.get(t - j), delta.get(t - j)));
			}
		}
		R min = lattice.zero();
		for (int i = j; i <= k; i++) {
			min = lattice.add(lattice.scalarMultiply(minParameters.get(i - j), basis.get(i)), min);
		}
		return min;
	}

//		List<S> orthogonal = gramSchmidt(embedList(basis, lattice));
//		List<List<T>> coefficients = new ArrayList<>();
//		for (int i = 0; i < basis.size(); i++) {
//			coefficients.add(new ArrayList<>());
//			for (int j = 0; j < basis.size(); j++) {
//				coefficients.get(i).add(f.zero());
//			}
//		}
//		computeLLLCoefficients(basis, orthogonal, coefficients, lattice);
//		int k = 1;
//		while (k < basis.size()) {
//			boolean adjusted = false;
//			for (int j = k - 1; j >= 0; j--) {
//				if (!r.close(f.value(coefficients.get(k).get(j)),half) && f.value(coefficients.get(k).get(j)).compareTo(half) > 0) {
//					adjusted = true;
//					basis.set(k, lattice.subtract(basis.get(k),
//							lattice.scalarMultiply(round(coefficients.get(k).get(j)), basis.get(j))));
//					orthogonal = gramSchmidt(embedList(basis, lattice));
//					computeLLLCoefficients(basis, orthogonal, coefficients, lattice);
//				}
//			}
//			Real c = f.value(coefficients.get(k).get(k - 1));
//			Real ksquare = f.value(innerProduct(orthogonal.get(k), orthogonal.get(k)));
//			Real km1square = f.value(innerProduct(orthogonal.get(k - 1), orthogonal.get(k - 1)));
//			Real rhs = r.multiply(r.subtract(delta, r.multiply(c, c)), km1square);
//			if (r.close(ksquare, rhs) || ksquare.compareTo(rhs) >= 0) {
//				k = adjusted ? k : k+1;
//			} else {
//				R tmp = basis.get(k);
//				basis.set(k, basis.get(k - 1));
//				basis.set(k - 1, tmp);
//				orthogonal = gramSchmidt(embedList(basis, lattice));
//				computeLLLCoefficients(basis, orthogonal, coefficients, lattice);
//				k = Math.max(k - 1, 1);
//			}
//		}
//		return basis;

//	private <R extends Element<R>> List<S> embedList(List<R> list, Lattice<R, T, S> lattice) {
//		List<S> result = new ArrayList<>();
//		for (R t : list) {
//			result.add(lattice.embedding(t));
//		}
//		return result;
//	}

//	private <R extends Element<R>> void computeLLLCoefficients(List<R> basis, List<S> orthogonal,
//			List<List<T>> coefficients, Lattice<R, T, S> lattice) {
//		ValueField<T> f = getValueField();
//		for (int i = 0; i < basis.size(); i++) {
//			for (int j = 0; j < basis.size(); j++) {
//				coefficients.get(i).set(j, f.divide(innerProduct(lattice.embedding(basis.get(i)), orthogonal.get(j)),
//						innerProduct(orthogonal.get(j), orthogonal.get(j))));
//			}
//		}
//	}

	@Override
	public <R extends Element<R>> R closestLatticePoint(S t, Lattice<R, T, S> lattice) {
		if (lattice.getModuleGenerators().isEmpty()) {
			return lattice.zero();
		}
		Matrix<T> asMatrix = lattice.generatorsAsMatrix();
		Matrix<T> pseudoInverse = pseudoInverse(asMatrix);
		MatrixModule<T> matrixModule = pseudoInverse.getModule(getField());
		Vector<T> solved = matrixModule.multiply(pseudoInverse, asVector(t));
		List<IntE> integerSolved = new ArrayList<>();
		for (T coefficient : solved.asList()) {
			integerSolved.add(round(coefficient));
		}
		return lattice.fromVector(new Vector<>(integerSolved));
	}

	@Override
	public <R extends Element<R>> List<R> latticePointsInParallelotope(S edge, Lattice<R, T, S> lattice) {
		Vector<T> asVector = asVector(edge);
		ValueField<T> field = getValueField();
		DualVectorSpace<T, S> dualSpace = getDual();
		List<Dual<T, S>> constraints = new ArrayList<>();
		List<T> rhs = new ArrayList<>();
		for (int i = 0; i < dimension(); i++) {
			constraints.add(dualSpace.negative(dualSpace.getUnitVector(i + 1)));
			rhs.add(field.zero());
			constraints.add(dualSpace.getUnitVector(i + 1));
			rhs.add(asVector.get(i + 1));
		}
		return latticePointsInPolytope(new Polytope<>(this, constraints, rhs), lattice);
	}

	@Override
	public <R extends Element<R>> List<R> latticePointsInPolytope(Polytope<T, S> polytope, Lattice<R, T, S> lattice) {
		return latticePointsInPolytope(polytope, lattice, 0.75);
	}

	@Override
	public <R extends Element<R>> List<R> latticePointsInPolytope(Polytope<T, S> polytope, Lattice<R, T, S> lattice,
			double delta) {
		Integers z = Integers.z();
		DualVectorSpace<T, S> dual = getDual();
		List<Dual<T, S>> dualBasis = dual.getDualBasis(lattice.generatorsAsMatrix());
		IntE crossProductSize = z.one();
		List<IntE> minima = new ArrayList<>();
		List<IntE> maxima = new ArrayList<>();
		for (int i = 0; i < lattice.rank(); i++) {
			IntE minValue = asReal(polytope.minimize(dualBasis.get(i)).getValue()).round();
			IntE maxValue = asReal(polytope.maximize(dualBasis.get(i)).getValue()).round();
			crossProductSize = z.multiply(z.add(z.subtract(maxValue, minValue), z.one()), crossProductSize);
			minima.add(minValue);
			maxima.add(maxValue);
		}
		if (crossProductSize.compareTo(z.getInteger(65536)) > 0) {
			return latticePointsInPolytope(polytope, lattice, dualBasis, Collections.emptyList());
		}
		List<List<IntE>> ranges = new ArrayList<>();
		for (int i = 0; i < lattice.rank(); i++) {
			IntE minValue = minima.get(i);
			IntE maxValue = maxima.get(i);
			List<IntE> range = new ArrayList<>();
			for (IntE counter = minValue; counter.compareTo(maxValue) <= 0; counter = z.add(counter, z.one())) {
				range.add(counter);
			}
			ranges.add(range);
		}
		List<List<IntE>> possiblePoints = MiscAlgorithms.crossProduct(ranges);
		List<R> result = new ArrayList<>();
		for (List<IntE> possiblePoint : possiblePoints) {
			R point = lattice.fromVector(new Vector<>(possiblePoint));// <Fraction> asFractionPoint = new ArrayList<>();
			if (!polytope.contains(lattice.embedding(point))) {
				continue;
			}
			result.add(lattice.fromVector(new Vector<>(possiblePoint)));
		}
		return result;
	}

	private <R extends Element<R>> List<R> latticePointsInPolytope(Polytope<T, S> polytope, Lattice<R, T, S> lattice,
			List<Dual<T, S>> dualLatticeBasis, List<IntE> partialVector) {
		int coordinate = partialVector.size();
		if (coordinate == lattice.rank()) {
			R point = lattice.fromVector(new Vector<>(partialVector));
			if (polytope.contains(lattice.embedding(point))) {
				return Collections.singletonList(point);
			}
			return Collections.emptyList();
		}
		if (polytope.isEmpty()) {
			return Collections.emptyList();
		}
		Integers z = Integers.z();
		ValueField<T> field = getValueField();
		Reals r = field.getReals();
		List<R> result = new ArrayList<>();
		DualVectorSpace<T, S> dualSpace = getDual();
		Dual<T, S> dualVector = dualLatticeBasis.get(coordinate);
		Dual<T, S> negativeDualVector = dualSpace.negative(dualVector);
		T epsilon = field.exactness().equals(Exactness.EXACT) ? field.zero()
				: fromReal(r.getPowerOfTwo(-r.precision() / 4));
		IntE minValue = asReal(polytope.minimize(dualVector).getValue()).round();
		IntE maxValue = asReal(polytope.maximize(dualVector).getValue()).round();
		List<Dual<T, S>> constraints = new ArrayList<>();
		constraints.add(dualVector);
		constraints.add(negativeDualVector);
		List<IntE> nextPartialVector = new ArrayList<>();
		nextPartialVector.addAll(partialVector);
		nextPartialVector.add(z.zero());
		for (IntE value = minValue; value.compareTo(maxValue) <= 0; value = z.add(value, z.one())) {
			nextPartialVector.set(coordinate, value);
			List<T> values = new ArrayList<>();
			values.add(field.add(field.getInteger(value), epsilon));
			values.add(field.add(field.getInteger(z.negative(value)), epsilon));
			result.addAll(latticePointsInPolytope(polytope.addConstraints(constraints, values), lattice,
					dualLatticeBasis, nextPartialVector));
		}
		return result;
	}

	@Override
	public <R extends Element<R>> List<R> latticeVertexPointsInPolytope(Polytope<T, S> polytope,
			Lattice<R, T, S> lattice) {
		Set<R> latticeVertexPoints = new TreeSet<>();
		List<S> vertices = new ArrayList<>();
		boolean newPointsFound;
		Polytope<T, S> latticePolytope = polytope;
		while (true) {
			newPointsFound = false;
			for (Dual<T, S> constraint : latticePolytope.getConstraints()) {
				R extremalPoint = latticeLinearProgram(polytope, constraint, lattice);
				if (extremalPoint != null && !latticeVertexPoints.contains(extremalPoint)) {
					latticeVertexPoints.add(extremalPoint);
					newPointsFound = true;
					vertices.add(lattice.embedding(extremalPoint));
				}
			}
			if (!newPointsFound) {
				break;
			}
			latticePolytope = Polytope.fromVertices(this, vertices);
		}
		Set<R> result = new TreeSet<>();
		for (S vertex : latticePolytope.vertices()) {
			result.add(closestLatticePoint(vertex, lattice));
		}
		List<R> asList = new ArrayList<>();
		asList.addAll(result);
		return asList;
	}

	@Override
	public <R extends Element<R>> R latticeLinearProgram(Polytope<T, S> polytope, Dual<T, S> maximize,
			Lattice<R, T, S> lattice) {
		ValueField<T> field = getValueField();
		DualVectorSpace<T, S> dualSpace = getDual();
		Matrix<T> baseChange = lattice.generatorsAsMatrix();
		List<Dual<T, S>> dualBasis = dualSpace.getDualBasis(baseChange);
		Matrix<T> invertedBaseChange = matrixAlgebra().inverse(baseChange);
		Deque<Polytope<T, S>> problemQueue = new LinkedList<>();
		problemQueue.add(polytope);
		T bestValue = null;
		R bestSolution = null;
		T errorMargin;
		if (field.exactness().equals(Exactness.EXACT)) {
			errorMargin = field.zero();
		} else {
			errorMargin = field.power(field.getInteger(2), -96);
		}
		while (!problemQueue.isEmpty()) {
			Polytope<T, S> problem = problemQueue.poll();
			LinearProgramResult<T, S> relaxedSolution = problem.maximize(maximize);
			if (relaxedSolution.getStatus().equals(SimplexAlgorithmStatus.EMPTY_POLYTOPE)) {
				continue;
			}
			if (relaxedSolution.getStatus().equals(SimplexAlgorithmStatus.UNBOUNDED_SOLUTION)) {
				throw new ArithmeticException("Not a bounded polytope!");
			}
			if (!relaxedSolution.getStatus().equals(SimplexAlgorithmStatus.OKAY)) {
				System.err.println(relaxedSolution.getStatus());
				System.err.println(problem.getConstraints());
				System.err.println(problem.getRightHandSide());
				throw new ArithmeticException("Simplex failed!");
			}
			if (bestValue != null && bestValue.compareTo(relaxedSolution.getValue()) >= 0) {
				continue;
			}
			Vector<T> solutionInLatticeBase = matrixAlgebra().multiply(invertedBaseChange,
					asVector(relaxedSolution.getSolution()));
			if (isIntegerVector(solutionInLatticeBase)) {
				bestValue = relaxedSolution.getValue();
				bestSolution = closestLatticePoint(relaxedSolution.getSolution(), lattice);
				continue;
			}
			T maximalFractionalPart = field.zero();
			int indexOfMaximalFractionalPart = -1;
			for (int i = 0; i < dimension(); i++) {
				T coefficient = solutionInLatticeBase.get(i + 1);
				T fractionalPart = abs(field.subtract(coefficient, field.getInteger(round(coefficient))));
				if (fractionalPart.compareTo(maximalFractionalPart) > 0) {
					indexOfMaximalFractionalPart = i + 1;
					maximalFractionalPart = fractionalPart;
				}
			}
			List<Dual<T, S>> lowerConstraints = new ArrayList<>();
			List<T> lowerRhs = new ArrayList<>();
			lowerConstraints.add(dualBasis.get(indexOfMaximalFractionalPart - 1));
			lowerRhs.add(field.add(field.getInteger(roundDown(solutionInLatticeBase.get(indexOfMaximalFractionalPart))),
					errorMargin));
			problemQueue.add(problem.addConstraints(lowerConstraints, lowerRhs));
			List<Dual<T, S>> upperConstraints = new ArrayList<>();
			List<T> upperRhs = new ArrayList<>();
			upperConstraints.add(dualSpace.negative(dualBasis.get(indexOfMaximalFractionalPart - 1)));
			upperRhs.add(field.negative(field.subtract(
					field.getInteger(roundUp(solutionInLatticeBase.get(indexOfMaximalFractionalPart))), errorMargin)));
			problemQueue.add(problem.addConstraints(upperConstraints, upperRhs));
		}
		return bestSolution;
	}

	private IntE roundDown(T t) {
		ValueField<T> field = getValueField();
		IntE rounded = round(t);
		T embeded = field.getInteger(rounded);
		if (t.compareTo(embeded) >= 0) {
			return rounded;
		}
		Integers z = Integers.z();
		return z.subtract(rounded, z.one());
	}

	private IntE roundUp(T t) {
		ValueField<T> field = getValueField();
		IntE rounded = round(t);
		T embeded = field.getInteger(rounded);
		if (embeded.compareTo(t) >= 0) {
			return rounded;
		}
		Integers z = Integers.z();
		return z.add(rounded, z.one());

	}

	private T abs(T t) {
		if (t.compareTo(getValueField().zero()) >= 0) {
			return t;
		}
		return getValueField().negative(t);
	}

	private boolean isInteger(T t) {
		Reals r = Reals.r(128);
		Real threshold = r.getPowerOfTwo(-64);
		ValueField<T> field = getValueField();
		return r.abs(r.subtract(field.value(t), field.value(field.getInteger(round(t))))).compareTo(threshold) < 0;
	}

	private boolean isIntegerVector(Vector<T> t) {
		for (T c : t.asList()) {
			if (!isInteger(c)) {
				return false;
			}
		}
		return true;
	}

}
