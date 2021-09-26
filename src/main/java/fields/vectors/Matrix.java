package fields.vectors;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;

import fields.helper.AbstractElement;
import fields.interfaces.BilinearMap;
import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.interfaces.Module;
import fields.interfaces.Ring;

public class Matrix<T extends Element<T>> extends AbstractElement<Matrix<T>> {
	private T[][] elements;
	MatrixModule<T>.SmithNormalFormResult smithNormalFormResult = null;
	MatrixModule<T>.LDUPResult ldupResult = null;

	@SuppressWarnings("unchecked")
	public Matrix(List<List<T>> elements) {
		this.elements = (T[][]) Array.newInstance(elements.get(0).get(0).getClass(), elements.size(),
				elements.get(0).size());
		for (int i = 0; i < elements.size(); i++) {
			for (int j = 0; j < elements.get(i).size(); j++) {
				this.elements[i][j] = elements.get(i).get(j);
			}
		}
	}

	public static <T extends Element<T>> Matrix<T> fromRows(List<Vector<T>> rows) {
		List<List<T>> m = new ArrayList<>();
		for (int i = 0; i < rows.size(); i++) {
			m.add(rows.get(i).asList());
		}
		return new Matrix<T>(m);
	}

	public static <T extends Element<T>> Matrix<T> fromColumns(List<Vector<T>> columns) {
		List<List<T>> m = new ArrayList<>();
		for (int i = 0; i < columns.get(0).dimension(); i++) {
			m.add(new ArrayList<>());
			for (int j = 0; j < columns.size(); j++) {
				m.get(i).add(columns.get(j).get(i + 1));
			}
		}
		return new Matrix<T>(m);
	}

	public static <T extends Element<T>, S extends Element<S>, U extends Element<U>> Matrix<T> fromLinearMap(
			Module<T, S> domain, Module<T, U> codomain, MathMap<S, U> map) {
		List<Vector<T>> m = new ArrayList<>();
		for (S s : domain.getModuleGenerators()) {
			m.add(codomain.asVector(map.evaluate(s)));
		}
		return fromColumns(m);
	}

	public static <T extends Element<T>, S extends Element<S>> Matrix<T> fromEndomorphism(Module<T, S> domain,
			MathMap<S, S> map) {
		return fromLinearMap(domain, domain, map);
	}
	
	public static <T extends Element<T>, S extends Element<S>> Matrix<T> fromBilinearMap(Module<T, S> space, BilinearMap<S, T> map) {
		List<List<T>> elements = new ArrayList<>();
		for (S s1 : space.getModuleGenerators()) {
			List<T> row = new ArrayList<>();
			for (S s2 : space.getModuleGenerators()) {
				row.add(map.evaluate(s1, s2));
			}
			elements.add(row);
		}
		return new Matrix<>(elements);
	}

	public static <T extends Element<T>, S extends Element<S>> Matrix<S> mapMatrix(MathMap<T, S> map, Matrix<T> t) {
		List<List<S>> result = new ArrayList<>();
		for (int i = 0; i < t.rows(); i++) {
			List<S> rowResult = new ArrayList<>();
			for (int j = 0; j < t.columns(); j++) {
				rowResult.add(map.evaluate(t.entry(i + 1, j + 1)));
			}
			result.add(rowResult);
		}
		return new Matrix<>(result);
	}

	@SuppressWarnings("unchecked")
	public Matrix(T[][] elements) {
		this.elements = (T[][]) Array.newInstance(elements[0][0].getClass(), elements.length, elements[0].length);
		for (int i = 0; i < elements.length; i++) {
			for (int j = 0; j < elements[i].length; j++) {
				this.elements[i][j] = elements[i][j];
			}
		}
	}

	public MatrixModule<T> getModule(Ring<T> base) {
		return new MatrixModule<>(base, rows(), columns());
	}

	public boolean isSquare() {
		return rows() == columns();
	}

	public int rows() {
		return elements.length;
	}

	public int columns() {
		return elements[0].length;
	}

	public Vector<T> row(int row) {
		List<T> r = new ArrayList<>();
		for (int i = 0; i < columns(); i++) {
			r.add(elements[row - 1][i]);
		}
		return new Vector<T>(r);
	}

	public Vector<T> column(int column) {
		List<T> r = new ArrayList<>();
		for (int i = 0; i < rows(); i++) {
			r.add(elements[i][column - 1]);
		}
		return new Vector<T>(r);
	}

	public T entry(int row, int column) {
		return elements[row - 1][column - 1];
	}

	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		int max = -1;
		List<List<String>> strings = new ArrayList<>();
		for (int i = 0; i < rows(); i++) {
			strings.add(new ArrayList<String>());
			for (int j = 0; j < columns(); j++) {
				String here = elements[i][j].toString();
				strings.get(i).add(here);
				if (max < here.length()) {
					max = here.length();
				}
			}
		}

		for (int i = 0; i < rows(); i++) {
			buf.append("[");
			for (int j = 0; j < columns(); j++) {
				if (j != 0) {
					buf.append("  ");
				}
				String here = strings.get(i).get(j);
				buf.append(" ".repeat(max - here.length()));
				buf.append(here);
			}
			buf.append("]\n");
		}
		return buf.toString();
	}

	@Override
	public int compareTo(Matrix<T> o) {
		if (o.rows() != this.rows()) {
			return o.rows() - this.rows();
		}
		if (o.columns() != this.columns()) {
			return o.columns() - this.columns();
		}
		for (int i = 0; i < elements.length; i++) {
			for (int j = 0; j < elements[i].length; j++) {
				int cmp = this.elements[i][j].compareTo(o.elements[i][j]);
				if (cmp != 0) {
					return cmp;
				}
			}
		}
		return 0;
	}

}
