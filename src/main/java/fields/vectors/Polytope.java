package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.InnerProductSpace.OrthogonalSimilarResult;
import fields.interfaces.MathSet;
import fields.interfaces.RealInnerProductSpace;
import fields.interfaces.RealInnerProductSpace.LinearProgramResult;
import fields.interfaces.RealInnerProductSpace.SimplexAlgorithmStatus;
import fields.interfaces.ValueField;
import fields.vectors.DualVectorSpace.Dual;

public class Polytope<T extends Element<T>, S extends Element<S>> implements MathSet<S> {
	private RealInnerProductSpace<T, S> space;
	private List<Dual<T, S>> constraints;
	private List<T> rhs;
	private List<S> vertices;
	private Tableau knownVertex;
	private boolean empty;
	private boolean deficient;

	public Polytope(RealInnerProductSpace<T, S> space, List<Dual<T, S>> constraints, List<T> rhs) {
		this.space = space;
		this.constraints = constraints;
		this.rhs = rhs;
		if (constraints.size() != rhs.size()) {
			throw new ArithmeticException("Wrong dimensions!");
		}
	}

	private static class Face<T extends Element<T>, S extends Element<S>> extends AbstractElement<Face<T, S>> {
		private SortedSet<S> vertices;
		private Dual<T, S> constraint;
		private T value;

		private Face(SortedSet<S> vertices, Dual<T, S> constraint, T value) {
			this.vertices = vertices;
			this.constraint = constraint;
			this.value = value;
		}

		@Override
		public int compareTo(Face<T, S> o) {
			Iterator<S> it = vertices.iterator();
			Iterator<S> oIt = o.vertices.iterator();
			while (it.hasNext()) {
				S next = it.next();
				S oNext = oIt.next();
				int cmp = next.compareTo(oNext);
				if (cmp != 0) {
					return cmp;
				}
			}
			return 0;
		}

		public String toString() {
			return vertices.toString() + " " + constraint + " " + value;
		}

	}

	private static <T extends Element<T>, S extends Element<S>> List<Face<T, S>> getFacesOfSimplex(
			RealInnerProductSpace<T, S> space, List<S> vertices) {
		DualVectorSpace<T, S> dual = space.getDual();
		List<Face<T, S>> result = new ArrayList<>();
		for (int i = 0; i < vertices.size(); i++) {
			S leftOut = vertices.get(i);
			SortedSet<S> faceVertices = new TreeSet<>();
			faceVertices.addAll(vertices);
			faceVertices.remove(leftOut);
			S baseVector = faceVertices.first();
			List<S> rebased = new ArrayList<>();
			for (S vertex : faceVertices) {
				if (!vertex.equals(baseVector)) {
					rebased.add(space.subtract(vertex, baseVector));
				}
			}
			SubRealInnerProductSpace<T, S> subSpace = new SubRealInnerProductSpace<>(space, rebased, false);
			SubRealInnerProductSpace<T, S> orthogonal = subSpace.orthogonalComplement();
			if (orthogonal.dimension() != 1) {
				throw new ArithmeticException("Not a full dimensional simplex");
			}
			Dual<T, S> constraint = dual.canonicalIsomorphism(orthogonal.getBasis().get(0));
			T value = constraint.evaluate(baseVector);
			if (constraint.evaluate(leftOut).compareTo(value) > 0) {
				constraint = dual.negative(constraint);
				value = space.getValueField().negative(value);
			}
			result.add(new Face<>(faceVertices, constraint, value));
		}
		return result;
	}

	private static <T extends Element<T>, S extends Element<S>> T getDistanceToSimplex(
			RealInnerProductSpace<T, S> space, List<S> vertices, S newPoint) {
		S basePoint = vertices.get(0);
		List<S> rebased = new ArrayList<>();
		for (S vertex : vertices) {
			if (vertex.equals(basePoint)) {
				continue;
			}
			rebased.add(space.subtract(vertex, basePoint));
		}
		rebased = space.gramSchmidt(rebased);
		S orthogonalNewPoint = newPoint;
		for (S orthogonal : rebased) {
			orthogonalNewPoint = space.subtract(orthogonalNewPoint,
					space.scalarMultiply(space.getField().divide(space.innerProduct(orthogonal, orthogonalNewPoint),
							space.innerProduct(orthogonal, orthogonal)), orthogonal));
		}
		return space.innerProduct(orthogonalNewPoint, orthogonalNewPoint);
	}

	public static <T extends Element<T>, S extends Element<S>> Polytope<T, S> fromVertices(
			RealInnerProductSpace<T, S> space, List<S> vertices) {
		DualVectorSpace<T, S> dual = space.getDual();
		ValueField<T> field = space.getValueField();
		if (vertices.isEmpty()) {
			List<Dual<T, S>> constraints = new ArrayList<>();
			List<T> values = new ArrayList<>();
			constraints.add(dual.getUnitVector(1));
			constraints.add(dual.negative(dual.getUnitVector(1)));
			values.add(field.zero());
			values.add(field.one());
			return new Polytope<>(space, constraints, values);
		}
		S baseVector = vertices.get(0);
		List<S> rebase = new ArrayList<>();
		for (int i = 1; i < vertices.size(); i++) {
			rebase.add(space.subtract(vertices.get(i), baseVector));
		}
		SubRealInnerProductSpace<T, S> subVectorSpace = new SubRealInnerProductSpace<>(space, rebase);
		if (subVectorSpace.dimension() < space.dimension()) {
			Polytope<T, S> polytope = fromVertices(subVectorSpace, vertices);
			SubRealInnerProductSpace<T, S> complement = subVectorSpace.orthogonalComplement();
			List<Dual<T, S>> constraints = new ArrayList<>();
			List<T> values = new ArrayList<>();
			for (Dual<T, S> constraint : polytope.getConstraints()) {
				constraints.add(dual.canonicalIsomorphism(constraint.dual()));
			}
			values.addAll(polytope.getRightHandSide());
			List<S> basis = complement.gramSchmidt(complement.getBasis());
			for (S basisVector : basis) {
				Dual<T, S> constraint = dual.canonicalIsomorphism(basisVector);
				T value = constraint.evaluate(baseVector);
				constraints.add(constraint);
				values.add(value);
				constraints.add(dual.negative(constraint));
				values.add(field.negative(value));
			}
			return new Polytope<>(space, constraints, values);
		}
		List<S> simplex = new ArrayList<>();
		S minLexVertex = null;
		S maxLexVertex = null;
		for (S vertex : vertices) {
			if (minLexVertex == null || minLexVertex.compareTo(vertex) > 0) {
				minLexVertex = vertex;
			}
			if (maxLexVertex == null || maxLexVertex.compareTo(vertex) < 0) {
				maxLexVertex = vertex;
			}
		}
		simplex.add(minLexVertex);
		simplex.add(maxLexVertex);
		while (simplex.size() < space.dimension() + 1) {
			S maxDistanceVertex = null;
			T maxValue = field.zero();
			for (S vertex : vertices) {
				T value = getDistanceToSimplex(space, simplex, vertex);
				if (value.compareTo(maxValue) > 0) {
					maxDistanceVertex = vertex;
					maxValue = value;
				}
			}
			if (maxDistanceVertex == null) {
				throw new ArithmeticException("Not high enough dimension!");
			}
			simplex.add(maxDistanceVertex);
		}
		Deque<Face<T, S>> faceQueue = new LinkedList<>();
		faceQueue.addAll(getFacesOfSimplex(space, simplex));
		List<Face<T, S>> faceResult = new ArrayList<>();
		Real epsilon = field.value(field.power(field.getInteger(2), -64));
		while (!faceQueue.isEmpty()) {
			Face<T, S> face = faceQueue.poll();
			T maxValue = null;
			S maxVertex = null;
			for (S vertex : vertices) {
				T eval = face.constraint.evaluate(vertex);
				if (field.value(field.subtract(eval, face.value)).compareTo(epsilon) < 0) {
					continue;
				}
				if (eval.compareTo(face.value) > 0 && (maxValue == null || maxValue.compareTo(eval) < 0)) {
					maxValue = eval;
					maxVertex = vertex;
				}
			}
			if (maxValue == null) {
				faceResult.add(face);
				continue;
			}
			List<S> newSimplex = new ArrayList<>();
			newSimplex.addAll(face.vertices);
			newSimplex.add(maxVertex);
			List<Face<T, S>> faces = getFacesOfSimplex(space, newSimplex);
			for (Face<T, S> newFace : faces) {
				if (!face.equals(newFace) && !faceQueue.remove(newFace)) {
					faceQueue.add(newFace);
				}
			}
		}
		List<Dual<T, S>> constraints = new ArrayList<>();
		List<T> values = new ArrayList<>();
		Set<S> hull = new TreeSet<>();
		for (Face<T, S> face : faceResult) {
			constraints.add(face.constraint);
			values.add(face.value);
			hull.addAll(face.vertices);
		}
		Polytope<T, S> polytope = new Polytope<>(space, constraints, values);
		polytope.vertices = new ArrayList<>();
		polytope.vertices.addAll(hull);
		return polytope;
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		for (int i = 0; i < constraints.size(); i++) {
			builder.append(constraints.get(i));
			builder.append(".x <= ");
			builder.append(rhs.get(i));
			builder.append("\n");
		}
		return builder.toString();
	}

	public List<Dual<T, S>> getConstraints() {
		return constraints;
	}

	public List<T> getRightHandSide() {
		return rhs;
	}

	public boolean contains(S s) {
		for (int i = 0; i < constraints.size(); i++) {
			if (constraints.get(i).evaluate(s).compareTo(rhs.get(i)) > 0) {
				return false;
			}
		}
		return true;
	}

	public List<S> vertices() {
		if (vertices == null) {
			vertices = new ArrayList<>();
			if (isEmpty() || isDeficient()) {
				return vertices;
			}
			Set<S> uniqueVertices = new TreeSet<>();
			Set<Tableau> visited = new TreeSet<>();
			Deque<Tableau> queue = new LinkedList<>();
			queue.add(new Tableau(knownVertex));
			while (!queue.isEmpty()) {
				Tableau tableau = queue.poll();
				if (visited.contains(tableau)) {
					continue;
				}
				tableau.refreshTableau();
				S vertex = space.fromVector(tableau.vertexVector);
				visited.add(tableau);
				uniqueVertices.add(vertex);
				for (Tableau neighbor : neighbors(tableau)) {
					queue.add(neighbor);
				}
			}
			vertices.addAll(uniqueVertices);
		}
		return vertices;
	}

	public boolean isEmpty() {
		knownVertex();
		return empty;
	}

	public boolean isDeficient() {
		knownVertex();
		return deficient;
	}

	public S knownVertex() {
		if (knownVertex == null && !empty && !deficient) {
			ValueField<T> field = space.getValueField();
			List<Vector<T>> rows = new ArrayList<>();
			Set<Integer> independentIndeces = new TreeSet<>();
			for (Dual<T, S> constraint : constraints) {
				rows.add(space.asVector(constraint.dual()));
			}
			int dimension = space.dimension();
			Matrix<T> matrix = Matrix.fromRows(rows);
			Matrix<T> square = space.matrixAlgebra().multiply(space.conjugateTranspose(matrix), matrix);
			InnerProductSpace<T, S> asFiniteVectorSpace = space.withDimension(dimension);
			int rank = 0;
			Real epsilon = field.getReals().getPowerOfTwo(-field.getReals().precision() + 1);
			if (field.exactness().equals(Exactness.EXACT)) {
				rank = space.matrixAlgebra().rank(square);
			} else {
				OrthogonalSimilarResult<T> schurr = space.schurrForm(square);
				for (int i = 0; i < dimension; i++) {
					if (field.value(schurr.getOrthogonallySimilarMatrix().entry(i + 1, i + 1)).compareTo(epsilon) < 0) {
						break;
					}
					rank++;
				}
			}
			if (rank < dimension) {
				deficient = true;
				return null;
			}
			List<S> independentElements = new ArrayList<>();
			for (int i = 0; i < rows.size(); i++) {
				Vector<T> row = rows.get(i);
				S fromRow = asFiniteVectorSpace.fromVector(row);
				S reduced = fromRow;
				for (S basisElement : independentElements) {
					reduced = asFiniteVectorSpace.subtract(reduced,
							asFiniteVectorSpace.scalarMultiply(
									field.divide(asFiniteVectorSpace.innerProduct(basisElement, reduced),
											asFiniteVectorSpace.innerProduct(basisElement, basisElement)),
									basisElement));
				}
				if (asFiniteVectorSpace.valueNorm(reduced).compareTo(epsilon) < 0) {
					continue;
				}
				independentIndeces.add(i);
				independentElements.add(reduced);
				if (independentIndeces.size() == rank) {
					break;
				}
			}
			Tableau tableau = new Tableau(independentIndeces);
			SimplexPivot pivot;
			int counter = 0;
			while (true) {
				counter++;
				if (counter % (10 * space.dimension()) == 0) {
					tableau.refreshTableau();
				}
				pivot = tableau.selectRowPivot();
				if (!pivot.status.equals(SimplexAlgorithmStatus.OKAY)) {
					break;
				}
				tableau.simplexStep(pivot.pivotRow, pivot.pivotCol);
				tableau.perturbOptimize();
			}
			if (pivot.status.equals(SimplexAlgorithmStatus.UNBOUNDED_SOLUTION)) {
				empty = true;
				return null;
			} else if (!pivot.status.equals(SimplexAlgorithmStatus.OPTIMAL)) {
				throw new ArithmeticException("Did not find a vertex in a non-empty polytope!");
			}
			Set<Integer> indeces = new TreeSet<>();
			indeces.addAll(tableau.vertexIndeces.values());
			knownVertex = new Tableau(indeces);
		}
		return space.fromVector(knownVertex.vertexVector);
	}

	private List<Tableau> neighbors(Tableau tableau) {
		ValueField<T> field = space.getValueField();
		List<Tableau> result = new ArrayList<>();
		T epsilon = field.exactness().equals(Exactness.EXACT) ? field.zero()
				: field.negative(field.power(field.getInteger(2), -field.getReals().precision() + 1));
		for (int j = 0; j < space.dimension(); j++) {
			int pivotRow = -1;
			T minCoeff = field.zero();
			for (int i = 0; i < tableau.tableauMatrix.rows(); i++) {
				if (tableau.tableauRhs.get(i + 1).compareTo(field.zero()) > 0) {
					return Collections.emptyList();
				}
				if (tableau.tableauMatrix.entry(i + 1, j + 1).compareTo(epsilon) >= 0) {
					continue;
				}
				T charCoeff = field.divide(tableau.tableauRhs.get(i + 1), tableau.tableauMatrix.entry(i + 1, j + 1));
				if (pivotRow < 0 || charCoeff.compareTo(minCoeff) < 0) {
					pivotRow = i + 1;
					minCoeff = charCoeff;
				}
			}
			if (pivotRow < 0) {
				continue;
			}
			Tableau copy = new Tableau(tableau);
			copy.simplexStep(pivotRow, j + 1);
			result.add(copy);
		}
		return result;
	}

	public S vertex(Set<Integer> indeces) {
		List<Vector<T>> rows = new ArrayList<>();
		List<T> restrictedRhs = new ArrayList<>();
		for (int index : indeces) {
			rows.add(space.asVector(constraints.get(index).dual()));
			restrictedRhs.add(rhs.get(index));
		}
		Matrix<T> matrix = Matrix.fromRows(rows);
		return space.fromVector(space.matrixAlgebra().solve(matrix, new Vector<>(restrictedRhs)));
	}

	public LinearProgramResult<T, S> maximize(Dual<T, S> linearForm) {
		if (isEmpty()) {
			return new LinearProgramResult<>(SimplexAlgorithmStatus.EMPTY_POLYTOPE);
		}
		if (isDeficient()) {
			List<S> duals = new ArrayList<>();
			for (Dual<T, S> constraint : constraints) {
				duals.add(constraint.dual());
			}
			SubRealInnerProductSpace<T, S> subVectorSpace = new SubRealInnerProductSpace<>(space, duals);
			SubRealInnerProductSpace<T, S> orthogonalSpace = subVectorSpace.orthogonalComplement();
			LinearProgramResult<T, S> restrictedProgram = new Polytope<>(subVectorSpace, constraints, rhs)
					.maximize(space.getDual().canonicalIsomorphism(subVectorSpace.project(linearForm.dual())));
			if (restrictedProgram.getStatus().equals(SimplexAlgorithmStatus.EMPTY_POLYTOPE)) {
				return new LinearProgramResult<>(SimplexAlgorithmStatus.EMPTY_POLYTOPE);
			}
			for (S orthogonalBasisVector : orthogonalSpace.getBasis()) {
				if (!linearForm.evaluate(orthogonalBasisVector).equals(space.getField().zero())) {
					return new LinearProgramResult<>(SimplexAlgorithmStatus.UNBOUNDED_SOLUTION);
				}
			}
			return restrictedProgram;

		}
		return maximize(new Tableau(knownVertex, linearForm));
	}

	public LinearProgramResult<T, S> minimize(Dual<T, S> linearForm) {
		LinearProgramResult<T, S> negativeResult = maximize(space.getDual().negative(linearForm));
		if (!negativeResult.getStatus().equals(SimplexAlgorithmStatus.OKAY)) {
			return negativeResult;
		}
		return new LinearProgramResult<>(negativeResult.getSolution(),
				space.getField().negative(negativeResult.getValue()), negativeResult.getVertexIndeces());
	}

//	private LinearProgramResult<T, S> maximize(Dual<T, S> linearForm, Set<Integer> knownVertex) {
//		Tableau tableau;
//		try {
//			tableau = new Tableau(knownVertex, linearForm);
//		} catch (ArithmeticException e) {
//			return new LinearProgramResult<>(SimplexAlgorithmStatus.NOT_FEASIBLE_START);
//		}
//		return maximize(tableau);
//	}

	private LinearProgramResult<T, S> maximize(Tableau tableau) {
		int counter = 0;
		while (true) {
			counter++;
			if (counter % (10 * space.dimension()) == 0) {
				tableau.refreshTableau();
			}
			SimplexPivot pivot = tableau.selectSimplexPivot();
			if (pivot.status == SimplexAlgorithmStatus.OPTIMAL) {
				tableau.refreshTableau();
				Set<Integer> vertexIndecesSet = new TreeSet<>();
				vertexIndecesSet.addAll(tableau.vertexIndeces.values());
				return new LinearProgramResult<>(space.fromVector(tableau.vertexVector), tableau.tableauValue,
						vertexIndecesSet);
			} else if (pivot.status == SimplexAlgorithmStatus.UNBOUNDED_SOLUTION
					|| pivot.status == SimplexAlgorithmStatus.NOT_FEASIBLE_START) {
				return new LinearProgramResult<>(pivot.status);
			}
			tableau.simplexStep(pivot.pivotRow, pivot.pivotCol);
		}
	}

	private class SimplexPivot {
		private int pivotRow;
		private int pivotCol;
		private SimplexAlgorithmStatus status;

		private SimplexPivot(SimplexAlgorithmStatus status) {
			this.status = status;
		}

		private SimplexPivot(int pivotRow, int pivotCol) {
			this.status = SimplexAlgorithmStatus.OKAY;
			this.pivotRow = pivotRow;
			this.pivotCol = pivotCol;
		}
	}

	private class Tableau implements Comparable<Tableau> {
		private Matrix<T> tableauMatrix;
		private Matrix<T> tableauOptimize;
		private Vector<T> tableauRhs;
		private T tableauValue;
		private Matrix<T> invertedVertexMatrix;
		private Vector<T> vertexVector;
		private Map<Integer, Integer> vertexIndeces;
		private Map<Integer, Integer> nonVertexIndeces;
		private Dual<T, S> linearForm;

		private Tableau(Set<Integer> knownVertex, Dual<T, S> linearForm) {
			this.linearForm = linearForm;
			init(knownVertex);
		}

		private Tableau(Set<Integer> knownVertex) {
			int dimension = space.dimension();
			List<T> optimize = new ArrayList<>();
			ValueField<T> field = space.getValueField();
			for (int i = 0; i < dimension; i++) {
				optimize.add(field.one());
			}
			int equations = constraints.size();
			List<Vector<T>> vertexRows = new ArrayList<>();
			for (int i = 0; i < equations; i++) {
				Dual<T, S> constraint = constraints.get(i);
				if (knownVertex.contains(i)) {
					vertexRows.add(space.asVector(constraint.dual()));
				}
			}
			Matrix<T> vertexMatrix = Matrix.fromRows(vertexRows);
			Vector<T> optimization = new Vector<>(optimize);
			Vector<T> inBase = space.matrixAlgebra().multiply(space.conjugateTranspose(vertexMatrix), optimization);
			this.linearForm = space.getDual().canonicalIsomorphism(space.fromVector(inBase));
			init(knownVertex);
//			int dimension = space.dimension();
//			int equations = constraints.size();
//			List<Vector<T>> vertexRows = new ArrayList<>();
//			List<Vector<T>> nonVertexRows = new ArrayList<>();
//			List<T> vertexRhsList = new ArrayList<>();
//			List<T> nonVertexRhsList = new ArrayList<>();
//			vertexIndeces = new TreeMap<>();
//			nonVertexIndeces = new TreeMap<>();
//			for (int i = 0; i < equations; i++) {
//				Dual<T, S> constraint = constraints.get(i);
//				if (knownVertex.contains(i)) {
//					vertexIndeces.put(vertexRows.size(), i);
//					vertexRows.add(space.asVector(constraint.dual()));
//					vertexRhsList.add(rhs.get(i));
//				} else {
//					nonVertexIndeces.put(nonVertexRows.size(), i);
//					nonVertexRows.add(space.asVector(constraint.dual()));
//					nonVertexRhsList.add(rhs.get(i));
//				}
//			}
//			Matrix<T> vertexMatrix = Matrix.fromRows(vertexRows);
//			Vector<T> vertexRhs = new Vector<>(vertexRhsList);
//			if (!space.matrixAlgebra().isUnit(vertexMatrix)) {
//				throw new ArithmeticException("Wrong dimension!");
//			}
//			invertedVertexMatrix = space.matrixAlgebra().inverse(vertexMatrix);
//			vertexVector = space.matrixAlgebra().multiply(invertedVertexMatrix, vertexRhs);
//			// S vertex = space.fromVector(vertexVector);
//			MatrixModule<T> nonVertexModule = new MatrixModule<>(space.getField(), equations - dimension, dimension);
//			FiniteVectorSpace<T> nonVertexVectorSpace = new FiniteVectorSpace<>(space.getField(),
//					equations - dimension);
//			Matrix<T> nonVertexMatrix = Matrix.fromRows(nonVertexRows);
//			tableauMatrix = nonVertexModule.multiply(nonVertexMatrix, invertedVertexMatrix);
//			List<T> optimize = new ArrayList<>();
//			ValueField<T> field = space.getValueField();
//			for (int i = 0; i < dimension; i++) {
//				optimize.add(field.one());
//			}
//			tableauOptimize = new Matrix<>(Collections.singletonList(optimize));
//			tableauRhs = nonVertexVectorSpace.subtract(nonVertexModule.multiply(tableauMatrix, vertexRhs),
//					new Vector<>(nonVertexRhsList));
//			List<T> roundedRhs = new ArrayList<>();
//			if (!space.getValueField().exactness().equals(Exactness.EXACT)) {
//				Reals r = field.getReals();
//				Real epsilon = r.getPowerOfTwo(-r.precision() + dimension);
//				for (T rhs : tableauRhs.asList()) {
//					if (field.value(rhs).compareTo(epsilon) < 0) {
//						roundedRhs.add(field.zero());
//					} else {
//						roundedRhs.add(rhs);
//					}
//				}
//				tableauRhs = new Vector<>(roundedRhs);
//			}
//			MatrixModule<T> optimizeModule = new MatrixModule<>(space.getField(), 1, dimension);
//			tableauValue = optimizeModule.multiply(optimizeModule.multiply(tableauOptimize, vertexMatrix), vertexVector)
//					.get(1);
		}

		private Tableau(Tableau tableau, Dual<T, S> linearForm) {
			this(tableau);
			this.linearForm = linearForm;
			this.tableauOptimize = new MatrixModule<>(space.getField(), 1, space.dimension())
					.multiply(linearForm.asMatrix(), tableau.invertedVertexMatrix);
			this.tableauValue = linearForm.evaluate(space.fromVector(vertexVector));
		}

		private Tableau(Tableau tableau) {
			this.linearForm = tableau.linearForm;
			this.tableauMatrix = tableau.tableauMatrix;
			this.tableauOptimize = tableau.tableauOptimize;
			this.tableauRhs = tableau.tableauRhs;
			this.tableauValue = tableau.tableauValue;
			this.invertedVertexMatrix = tableau.invertedVertexMatrix;
			this.vertexVector = tableau.vertexVector;
			this.vertexIndeces = new TreeMap<>();
			this.vertexIndeces.putAll(tableau.vertexIndeces);
			this.nonVertexIndeces = new TreeMap<>();
			this.nonVertexIndeces.putAll(tableau.nonVertexIndeces);
		}

		private void init(Set<Integer> knownVertex) {
			int dimension = space.dimension();
			int equations = constraints.size();
			List<Vector<T>> vertexRows = new ArrayList<>();
			List<Vector<T>> nonVertexRows = new ArrayList<>();
			List<T> vertexRhsList = new ArrayList<>();
			List<T> nonVertexRhsList = new ArrayList<>();
			vertexIndeces = new TreeMap<>();
			nonVertexIndeces = new TreeMap<>();
			for (int i = 0; i < equations; i++) {
				Dual<T, S> constraint = constraints.get(i);
				if (knownVertex.contains(i)) {
					vertexIndeces.put(vertexRows.size(), i);
					vertexRows.add(space.asVector(constraint.dual()));
					vertexRhsList.add(rhs.get(i));
				} else {
					nonVertexIndeces.put(nonVertexRows.size(), i);
					nonVertexRows.add(space.asVector(constraint.dual()));
					nonVertexRhsList.add(rhs.get(i));
				}
			}
			Matrix<T> vertexMatrix = Matrix.fromRows(vertexRows);
			Vector<T> vertexRhs = new Vector<>(vertexRhsList);
			if (!space.matrixAlgebra().isUnit(vertexMatrix)) {
				throw new ArithmeticException("Wrong dimension!");
			}
			invertedVertexMatrix = space.matrixAlgebra().inverse(vertexMatrix);
			vertexVector = space.matrixAlgebra().multiply(invertedVertexMatrix, vertexRhs);
			S vertex = space.fromVector(vertexVector);
			MatrixModule<T> nonVertexModule = new MatrixModule<>(space.getField(), equations - dimension, dimension);
			FiniteVectorSpace<T> nonVertexVectorSpace = new FiniteVectorSpace<>(space.getField(),
					equations - dimension);
			Matrix<T> nonVertexMatrix = Matrix.fromRows(nonVertexRows);
			tableauMatrix = nonVertexModule.multiply(nonVertexMatrix, invertedVertexMatrix);
			tableauOptimize = new MatrixModule<>(space.getField(), 1, dimension).multiply(linearForm.asMatrix(),
					invertedVertexMatrix);
			tableauValue = linearForm.evaluate(vertex);
			tableauRhs = nonVertexVectorSpace.subtract(nonVertexModule.multiply(tableauMatrix, vertexRhs),
					new Vector<>(nonVertexRhsList));
			perturbRhs();
		}

		private void refreshTableau() {
			if (!space.getField().exactness().equals(Exactness.EXACT)) {
				Set<Integer> knownVertex = new TreeSet<>();
				knownVertex.addAll(vertexIndeces.values());
				init(knownVertex);
			}
		}

		@SuppressWarnings("unchecked")
		@Override
		public boolean equals(Object o) {
			Tableau other = (Tableau) o;
			return vertexIndeces.equals(other.vertexIndeces);
		}

		public int compareTo(Tableau other) {
			if (vertexIndeces.size() != other.vertexIndeces.size()) {
				throw new ArithmeticException("Comparing wrong sized tableaus!");
			}
			Set<Integer> vertexIndeces = new TreeSet<>();
			vertexIndeces.addAll(this.vertexIndeces.values());
			Set<Integer> otherVertexIndeces = new TreeSet<>();
			otherVertexIndeces.addAll(other.vertexIndeces.values());
			Iterator<Integer> it = vertexIndeces.iterator();
			Iterator<Integer> otherIt = otherVertexIndeces.iterator();
			while (it.hasNext()) {
				int cmp = it.next().compareTo(otherIt.next());
				if (cmp != 0) {
					return cmp;
				}
			}
			return 0;
		}

		private SimplexPivot selectSimplexPivot() {
			int pivotCol = -1;
			int pivotRow = -1;
			ValueField<T> field = space.getValueField();
			T epsilon = field.negative(field.exactness().equals(Exactness.EXACT) ? field.zero()
					: field.power(field.getInteger(2), -field.getReals().precision() + 1));
			T minCoeff = field.zero();
			for (int j = 0; j < space.dimension(); j++) {
				if (tableauOptimize.entry(1, j + 1).compareTo(epsilon) < 0) {
					pivotCol = j + 1;
					break;
				}
			}
			if (pivotCol < 0) {
				return new SimplexPivot(SimplexAlgorithmStatus.OPTIMAL);
			}
			for (int i = 0; i < tableauMatrix.rows(); i++) {
				if (tableauRhs.get(i + 1).compareTo(field.zero()) > 0) {
					System.err.println(vertexIndeces);
					System.err.println(nonVertexIndeces);
					System.err.println(tableauOptimize);
					System.err.println(tableauMatrix);
					System.err.println(tableauRhs);
					return new SimplexPivot(SimplexAlgorithmStatus.NOT_FEASIBLE_START);
				}
				if (tableauMatrix.entry(i + 1, pivotCol).compareTo(epsilon) >= 0) {
					continue;
				}
				T charCoeff = field.divide(tableauRhs.get(i + 1), tableauMatrix.entry(i + 1, pivotCol));
				if (pivotRow < 0 || charCoeff.compareTo(minCoeff) < 0) {
					pivotRow = i + 1;
					minCoeff = charCoeff;
				}
			}
			if (pivotRow < 0) {
				return new SimplexPivot(SimplexAlgorithmStatus.UNBOUNDED_SOLUTION);
			}
			return new SimplexPivot(pivotRow, pivotCol);
		}

		private SimplexPivot selectRowPivot() {
			int pivotCol = -1;
			int pivotRow = -1;
			ValueField<T> field = space.getValueField();
			T minCoeff = field.zero();
			for (int i = 0; i < tableauMatrix.rows(); i++) {
				if (tableauRhs.get(i + 1).compareTo(field.zero()) > 0) {
					pivotRow = i + 1;
					break;
				}
			}
			if (pivotRow < 0) {
				return new SimplexPivot(SimplexAlgorithmStatus.OPTIMAL);
			}
			T epsilon = field.exactness().equals(Exactness.EXACT) ? field.zero()
					: field.power(field.getInteger(2), -field.getReals().precision() + 1);

			for (int j = 0; j < tableauMatrix.columns(); j++) {
				if (tableauOptimize.entry(1, j + 1).compareTo(field.zero()) < 0) {
					System.err.println(vertexIndeces);
					System.err.println(nonVertexIndeces);
					System.err.println(tableauOptimize);
					System.err.println(tableauMatrix);
					System.err.println(tableauRhs);
					return new SimplexPivot(SimplexAlgorithmStatus.NOT_FEASIBLE_START);
				}
				if (tableauMatrix.entry(pivotRow, j + 1).compareTo(epsilon) <= 0) {
					continue;
				}
				T charCoeff = field.divide(tableauOptimize.entry(1, j + 1), tableauMatrix.entry(pivotRow, j + 1));
				if (pivotCol < 0 || charCoeff.compareTo(minCoeff) < 0) {
					pivotCol = j + 1;
					minCoeff = charCoeff;
				}
			}
			if (pivotCol < 0) {
				return new SimplexPivot(SimplexAlgorithmStatus.UNBOUNDED_SOLUTION);
			}
			return new SimplexPivot(pivotRow, pivotCol);
		}

		private void simplexStep(int pivotRow, int pivotCol) {
			int dimension = space.dimension();
			int nonVertexDimension = tableauRhs.dimension();
			Field<T> field = space.getField();
			T AMuNuInv = field.inverse(tableauMatrix.entry(pivotRow, pivotCol));
			List<List<T>> nextTableauMatrix = new ArrayList<>();
			for (int k = 0; k < nonVertexDimension; k++) {
				List<T> row = new ArrayList<>();
				for (int l = 0; l < dimension; l++) {
					if (k + 1 == pivotRow && l + 1 == pivotCol) {
						row.add(AMuNuInv);
					} else if (k + 1 == pivotRow) {
						row.add(field.negative(field.multiply(tableauMatrix.entry(pivotRow, l + 1), AMuNuInv)));
					} else if (l + 1 == pivotCol) {
						row.add(field.multiply(tableauMatrix.entry(k + 1, pivotCol), AMuNuInv));
					} else {
						row.add(field.subtract(tableauMatrix.entry(k + 1, l + 1), field.multiply(
								tableauMatrix.entry(k + 1, pivotCol), tableauMatrix.entry(pivotRow, l + 1), AMuNuInv)));
					}
				}
				nextTableauMatrix.add(row);
			}
			List<List<T>> nextVertexMatrix = new ArrayList<>();
			for (int k = 0; k < dimension; k++) {
				List<T> row = new ArrayList<>();
				for (int l = 0; l < dimension; l++) {
					if (l + 1 == pivotCol) {
						row.add(field.multiply(invertedVertexMatrix.entry(k + 1, pivotCol), AMuNuInv));
					} else {
						row.add(field.subtract(invertedVertexMatrix.entry(k + 1, l + 1),
								field.multiply(invertedVertexMatrix.entry(k + 1, pivotCol),
										tableauMatrix.entry(pivotRow, l + 1), AMuNuInv)));
					}
				}
				nextVertexMatrix.add(row);
			}
			List<T> nextTableauOptimize = new ArrayList<>();
			for (int l = 0; l < dimension; l++) {
				if (l + 1 == pivotCol) {
					nextTableauOptimize.add(field.multiply(tableauOptimize.entry(1, pivotCol), AMuNuInv));
				} else {
					nextTableauOptimize.add(field.subtract(tableauOptimize.entry(1, l + 1), field.multiply(
							tableauOptimize.entry(1, pivotCol), tableauMatrix.entry(pivotRow, l + 1), AMuNuInv)));
				}
			}
			List<T> nextTableauRhs = new ArrayList<>();
			for (int k = 0; k < nonVertexDimension; k++) {
				if (k + 1 == pivotRow) {
					nextTableauRhs.add(field.negative(field.multiply(tableauRhs.get(pivotRow), AMuNuInv)));
				} else {
					nextTableauRhs.add(field.subtract(tableauRhs.get(k + 1),
							field.multiply(tableauMatrix.entry(k + 1, pivotCol), tableauRhs.get(pivotRow), AMuNuInv)));
				}
			}
			List<T> nextVertex = new ArrayList<>();
			for (int k = 0; k < dimension; k++) {
				nextVertex.add(field.subtract(vertexVector.get(k + 1), field
						.multiply(invertedVertexMatrix.entry(k + 1, pivotCol), tableauRhs.get(pivotRow), AMuNuInv)));
			}
			T nextTableauValue = field.subtract(tableauValue,
					field.multiply(tableauOptimize.entry(1, pivotCol), tableauRhs.get(pivotRow), AMuNuInv));
			tableauMatrix = new Matrix<>(nextTableauMatrix);
			tableauOptimize = Matrix.fromRows(Collections.singletonList(new Vector<>(nextTableauOptimize)));
			tableauRhs = new Vector<>(nextTableauRhs);
			tableauValue = nextTableauValue;
			invertedVertexMatrix = new Matrix<>(nextVertexMatrix);
			vertexVector = new Vector<>(nextVertex);
			int nextNonVertexIndex = vertexIndeces.get(pivotCol - 1);
			int nextVertexIndex = nonVertexIndeces.get(pivotRow - 1);
			vertexIndeces.put(pivotCol - 1, nextVertexIndex);
			nonVertexIndeces.put(pivotRow - 1, nextNonVertexIndex);
			perturbRhs();
		}

		private void perturbRhs() {
			ValueField<T> field = space.getValueField();
			List<T> perturbedRhs = new ArrayList<>();
			if (!space.getValueField().exactness().equals(Exactness.EXACT)) {
				Reals r = field.getReals();
				Real epsilon = r.getPowerOfTwo(-r.precision() + 5);
				T pertubation = field.negative(field.power(field.getInteger(2), -r.precision() / 2));
				for (T rhs : tableauRhs.asList()) {
					if (field.value(rhs).compareTo(epsilon) < 0) {
						perturbedRhs.add(pertubation);
					} else {
						perturbedRhs.add(rhs);
					}
				}
				tableauRhs = new Vector<>(perturbedRhs);
			}
		}

		private void perturbOptimize() {
			ValueField<T> field = space.getValueField();
			List<T> perturbedOptimize = new ArrayList<>();
			if (!space.getValueField().exactness().equals(Exactness.EXACT)) {
				Reals r = field.getReals();
				Real epsilon = r.getPowerOfTwo(-r.precision() + 5);
				T pertubation = field.power(field.getInteger(2), -r.precision() / 2);
				for (T optimize : tableauOptimize.row(1).asList()) {
					if (field.value(optimize).compareTo(epsilon) < 0) {
						perturbedOptimize.add(pertubation);
					} else {
						perturbedOptimize.add(optimize);
					}
				}
				tableauOptimize = Matrix.fromRows(Collections.singletonList(new Vector<>(perturbedOptimize)));
			}
		}
	}

	@Override
	public Iterator<S> iterator() {
		throw new InfinityException();
	}

	@Override
	public Exactness exactness() {
		return space.exactness();
	}

	@Override
	public S getRandomElement() {
		List<S> vertices = vertices();
		S result = space.zero();
		Reals r = Reals.r(128);
		Real remaining = r.one();
		for (int i = 0; i < vertices.size() - 1; i++) {
			Real coefficient = r.getRandomElement(r.zero(), remaining);
			result = space.add(space.scalarMultiply(space.fromReal(coefficient), vertices.get(i)), result);
		}
		result = space.add(space.scalarMultiply(space.fromReal(remaining), vertices.get(vertices.size() - 1)), result);
		return result;
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

}
