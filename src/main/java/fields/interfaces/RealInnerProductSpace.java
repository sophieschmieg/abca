package fields.interfaces;

import java.util.List;
import java.util.Set;

import fields.floatingpoint.Reals.Real;
import fields.integers.Integers.IntE;
import fields.vectors.DualVectorSpace;
import fields.vectors.DualVectorSpace.Dual;
import fields.vectors.Polytope;

public interface RealInnerProductSpace<T extends Element<T>, S extends Element<S>> extends InnerProductSpace<T, S> {
	public T innerProduct(S s1, S s2);

	@Override
	public Real valueNorm(S s);

	public List<S> gramSchmidt(List<S> s);

	public List<S> normedGramSchmidt(List<S> s);

	public DualVectorSpace<T, S> getDual();

	public IntE round(T t);

	public T fromReal(Real r);

	public Real asReal(T t);

	public RealInnerProductSpace<T, S> withDimension(int dimension);

	public static enum SimplexAlgorithmStatus {
		OKAY, OPTIMAL, NOT_FEASIBLE_START, EMPTY_POLYTOPE, UNBOUNDED_SOLUTION;
	}

	public static class LinearProgramResult<T extends Element<T>, S extends Element<S>> {
		private SimplexAlgorithmStatus status;
		private S solution;
		private Set<Integer> vertexIndeces;
		private T value;

		public LinearProgramResult(SimplexAlgorithmStatus status) {
			this.status = status;
			this.solution = null;
			this.value = null;
		}

		public LinearProgramResult(S solution, T value, Set<Integer> vertexIndeces) {
			this.status = SimplexAlgorithmStatus.OKAY;
			this.solution = solution;
			this.value = value;
			this.vertexIndeces = vertexIndeces;
		}

		public SimplexAlgorithmStatus getStatus() {
			return status;
		}

		public S getSolution() {
			return solution;
		}

		public T getValue() {
			return value;
		}

		public Set<Integer> getVertexIndeces() {
			return vertexIndeces;
		}
	}

	public LinearProgramResult<T, S> linearProgram(Polytope<T, S> polytope, Dual<T, S> maximize);

	public <R extends Element<R>> List<R> latticeReduction(Lattice<R, T, S> lattice);

	public <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice);

	public <R extends Element<R>> List<R> latticeReduction(List<R> sublatticeBasis, Lattice<R, T, S> lattice, boolean isBasis);

	public <R extends Element<R>> R closestLatticePoint(S t, Lattice<R, T, S> lattice);

	public <R extends Element<R>> List<R> latticePointsInParallelotope(S edge, Lattice<R, T, S> lattice);

	public <R extends Element<R>> List<R> latticePointsInPolytope(Polytope<T, S> polytope, Lattice<R, T, S> lattice);

	public <R extends Element<R>> List<R> latticePointsInPolytope(Polytope<T, S> polytope, Lattice<R, T, S> lattice,
			double delta);

	public <R extends Element<R>> R latticeLinearProgram(Polytope<T, S> polytope, Dual<T, S> maximize,
			Lattice<R, T, S> lattice);

	public <R extends Element<R>> List<R> latticeVertexPointsInPolytope(Polytope<T, S> polytope,
			Lattice<R, T, S> lattice);

}
