package fields.vectors;

import java.util.ArrayList;
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
	public T conjugate(T t) {
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
	public <R extends Element<R>> List<R> latticeReduction(Lattice<R, T, S> lattice, double deltaAsDouble) {
		return latticeReduction(lattice.getModuleGenerators(), lattice, deltaAsDouble);
	}

	@Override
	public <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice) {
		return latticeReduction(sublatticeBasis, lattice, 0.75);
	}

	@Override
	public <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice,
			double deltaAsDouble) {
		ValueField<T> f = getValueField();
		Reals r = f.getReals();
		Real delta = r.getDouble(deltaAsDouble);
		Real half = r.getDouble(0.5);
		List<R> basis = new ArrayList<>();
		basis.addAll(sublatticeBasis);
		List<S> orthogonal = gramSchmidt(embedList(basis, lattice));
		List<List<T>> coefficients = new ArrayList<>();
		for (int i = 0; i < basis.size(); i++) {
			coefficients.add(new ArrayList<>());
			for (int j = 0; j < basis.size(); j++) {
				coefficients.get(i).add(f.zero());
			}
		}
		computeLLLCoefficients(basis, orthogonal, coefficients, lattice);
		int k = 1;
		while (k < basis.size()) {
			for (int j = k - 1; j >= 0; j--) {
				if (f.value(coefficients.get(k).get(j)).compareTo(half) > 0) {
					basis.set(k, lattice.subtract(basis.get(k),
							lattice.scalarMultiply(round(coefficients.get(k).get(j)), basis.get(j))));
					orthogonal = gramSchmidt(embedList(basis, lattice));
					computeLLLCoefficients(basis, orthogonal, coefficients, lattice);
				}
			}
			Real c = f.value(coefficients.get(k).get(k - 1));
			Real ksquare = f.value(innerProduct(orthogonal.get(k), orthogonal.get(k)));
			Real km1square = f.value(innerProduct(orthogonal.get(k - 1), orthogonal.get(k - 1)));
			Real rhs = r.multiply(r.subtract(delta, r.multiply(c, c)), km1square);
			if (r.close(ksquare, rhs) || ksquare.compareTo(rhs) >= 0) {
				k++;
			} else {
				R tmp = basis.get(k);
				basis.set(k, basis.get(k - 1));
				basis.set(k - 1, tmp);
				orthogonal = gramSchmidt(embedList(basis, lattice));
				computeLLLCoefficients(basis, orthogonal, coefficients, lattice);
				k = Math.max(k - 1, 1);
			}
		}
		return basis;
	}

	private <R extends Element<R>> List<S> embedList(List<R> list, Lattice<R, T, S> lattice) {
		List<S> result = new ArrayList<>();
		for (R t : list) {
			result.add(lattice.embedding(t));
		}
		return result;
	}

	private <R extends Element<R>> void computeLLLCoefficients(List<R> basis, List<S> orthogonal,
			List<List<T>> coefficients, Lattice<R, T, S> lattice) {
		ValueField<T> f = getValueField();
		for (int i = 0; i < basis.size(); i++) {
			for (int j = 0; j < basis.size(); j++) {
				coefficients.get(i).set(j, f.divide(innerProduct(lattice.embedding(basis.get(i)), orthogonal.get(j)),
						innerProduct(orthogonal.get(j), orthogonal.get(j))));
			}
		}
	}

	@Override
	public <R extends Element<R>> R closestLatticePoint(S t, Lattice<R, T, S> lattice) {
		return closestLatticePoint(t, lattice, 0.75);
	}

	@Override
	public <R extends Element<R>> R closestLatticePoint(S t, Lattice<R, T, S> lattice, double delta) {
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
		return latticePointsInParallelotope(edge, lattice, 0.75);
	}

	@Override
	public <R extends Element<R>> List<R> latticePointsInParallelotope(S edge, Lattice<R, T, S> lattice, double delta) {
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
		return latticePointsInPolytope(new Polytope<>(this, constraints, rhs), lattice, delta);
	}

	@Override
	public <R extends Element<R>> List<R> latticePointsInPolytope(Polytope<T, S> polytope, Lattice<R, T, S> lattice) {
		return latticePointsInPolytope(polytope, lattice, 0.75);
	}

	@Override
	public <R extends Element<R>> List<R> latticePointsInPolytope(Polytope<T, S> polytope, Lattice<R, T, S> lattice,
			double delta) {
		DualVectorSpace<T, S> dual = getDual();
		List<Dual<T, S>> dualBasis = dual.getDualBasis(lattice.generatorsAsMatrix());
		// List<R> vertices = latticeVertexPointsInPolytope(polytope, lattice, delta);
		Integers z = Integers.z();
//		Rationals q = Rationals.q();
//		FiniteRationalVectorSpace rationalSpace = new FiniteRationalVectorSpace(lattice.rank());
//		DualVectorSpace<Fraction, Vector<Fraction>> dualSpace = rationalSpace.getDual();
//		List<Vector<Fraction>> latticeVertices = new ArrayList<>();
//		for (R vertex : vertices) {
//			List<Fraction> rationalVector = new ArrayList<>();
//			Vector<IntE> integerVector = lattice.asVector(vertex);
//			for (IntE c : integerVector.asList()) {
//				rationalVector.add(q.getEmbedding(c));
//			}
//			latticeVertices.add(new Vector<>(rationalVector));
//		}
		// Polytope<Fraction, Vector<Fraction>> latticePolytope =
		// Polytope.fromVertices(rationalSpace, latticeVertices);
		List<List<IntE>> ranges = new ArrayList<>();
		for (int i = 0; i < lattice.rank(); i++) {
			// R minPoint = latticeLinearProgram(polytope, dual.negative(dualBasis.get(i)),
			// lattice);
			// R maxPoint = latticeLinearProgram(polytope, dualBasis.get(i), lattice);
			IntE minValue = asReal(polytope.minimize(dualBasis.get(i)).getValue()).round();// lattice.asVector(minPoint).get(i+1);
			IntE maxValue = asReal(polytope.maximize(dualBasis.get(i)).getValue()).round();// lattice.asVector(maxPoint).get(i+1);//latticePolytope.maximize(dualSpace.getUnitVector(i
																							// + 1)).getValue();
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
//			for (IntE c : possiblePoint) {
//				asFractionPoint.add(q.getInteger(c));
//			}
			if (!polytope.contains(lattice.embedding(point))) {
				continue;
			}
			result.add(lattice.fromVector(new Vector<>(possiblePoint)));
		}
		return result;
	}

	@Override
	public <R extends Element<R>> List<R> latticeVertexPointsInPolytope(Polytope<T, S> polytope,
			Lattice<R, T, S> lattice) {
		return latticeVertexPointsInPolytope(polytope, lattice, 0.75);
	}

	@Override
	public <R extends Element<R>> List<R> latticeVertexPointsInPolytope(Polytope<T, S> polytope,
			Lattice<R, T, S> lattice, double delta) {
		Set<R> latticeVertexPoints = new TreeSet<>();
		List<S> vertices = new ArrayList<>();
		boolean newPointsFound;
		Polytope<T, S> latticePolytope = polytope;
		while (true) {
			newPointsFound = false;
			for (Dual<T, S> constraint : latticePolytope.getConstraints()) {
				R extremalPoint = latticeLinearProgram(polytope, constraint, lattice, delta);
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
		return latticeLinearProgram(polytope, maximize, lattice, 0.75);
	}

	@Override
	public <R extends Element<R>> R latticeLinearProgram(Polytope<T, S> polytope, Dual<T, S> maximize,
			Lattice<R, T, S> lattice, double delta) {
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
			lowerConstraints.addAll(problem.getConstraints());
			List<T> lowerRhs = new ArrayList<>();
			lowerRhs.addAll(problem.getRightHandSide());
			lowerConstraints.add(dualBasis.get(indexOfMaximalFractionalPart - 1));
			lowerRhs.add(field.add(field.getInteger(roundDown(solutionInLatticeBase.get(indexOfMaximalFractionalPart))),
					errorMargin));
			problemQueue.add(new Polytope<>(this, lowerConstraints, lowerRhs));
			List<Dual<T, S>> upperConstraints = new ArrayList<>();
			upperConstraints.addAll(problem.getConstraints());
			List<T> upperRhs = new ArrayList<>();
			upperRhs.addAll(problem.getRightHandSide());
			upperConstraints.add(dualSpace.negative(dualBasis.get(indexOfMaximalFractionalPart - 1)));
			upperRhs.add(field.negative(field.subtract(
					field.getInteger(roundUp(solutionInLatticeBase.get(indexOfMaximalFractionalPart))), errorMargin)));
			problemQueue.add(new Polytope<>(this, upperConstraints, upperRhs));
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
