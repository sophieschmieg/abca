package fields.vectors;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.helper.AbstractModule;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Module;
import fields.interfaces.Ring;
import fields.vectors.TensorProductModule.TensorProduct;
import util.Pair;

public class TensorProductModule<T extends Element<T>, S extends Element<S>, U extends Element<U>>
		extends AbstractModule<T, TensorProduct<T, S, U>> {
	private static class TensorProductSummand<T extends Element<T>, S extends Element<S>, U extends Element<U>> {
		private T ringElement;
		private S firstGenerator;
		private U secondGenerator;
		
		private TensorProductSummand( T ringElement, S firstGenerator,
				U secondGenerator) {
			this.ringElement = ringElement;
			this.firstGenerator = firstGenerator;
			this.secondGenerator = secondGenerator;
		}

		public T getRingElement() {
			return ringElement;
		}

		public S getFirstGenerator() {
			return firstGenerator;
		}

		public U getSecondGenerator() {
			return secondGenerator;
		}
		
	}
	
	public static class TensorProduct<T extends Element<T>, S extends Element<S>, U extends Element<U>>
			extends AbstractElement<TensorProduct<T, S, U>> {
		private TensorProductModule<T, S, U> module;
		private List<TensorProductSummand<T, S, U>> summands;
	
		private TensorProduct(TensorProductModule<T, S, U> module, List<TensorProductSummand<T, S, U>> summands) {
			this.module = module;
			this.summands = summands;
		}
		
		
		
		@Override
		public int compareTo(TensorProduct<T, S, U> o) {
			// TODO Auto-generated method stub
			return 0;
		}

	}
	
	private Module<T, S> firstFactor;
	private Module<T, U> secondFactor;
	private Ring<T> base;
	
	public TensorProductModule(Module<T, S> firstFactor, Module<T, U> secondFactor) {
		if (!firstFactor.getRing().equals(secondFactor.getRing())) {
			throw new ArithmeticException("Different base rings!");
		}
		if (!firstFactor.isFree() || !secondFactor.isFree()) {
			throw new ArithmeticException("Factors are not free!");
		}
		this.base = firstFactor.getRing();
		this.firstFactor = firstFactor;
		this.secondFactor = secondFactor;
	}

	@Override
	public Ring<T> getRing() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TensorProduct<T, S, U> zero() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TensorProduct<T, S, U> add(TensorProduct<T, S, U> s1, TensorProduct<T, S, U> s2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TensorProduct<T, S, U> negative(TensorProduct<T, S, U> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TensorProduct<T, S, U> scalarMultiply(T t, TensorProduct<T, S, U> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isFree() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Ideal<T> annihilator() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isLinearIndependent(List<TensorProduct<T, S, U>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isGeneratingModule(List<TensorProduct<T, S, U>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<List<T>> nonTrivialCombinations(List<TensorProduct<T, S, U>> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<TensorProduct<T, S, U>> getModuleGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Vector<T> asVector(TensorProduct<T, S, U> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Exactness exactness() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TensorProduct<T, S, U> getRandomElement() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isFinite() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterator<TensorProduct<T, S, U>> iterator() {
		// TODO Auto-generated method stub
		return null;
	}
}
