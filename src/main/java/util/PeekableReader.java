package util;

import java.io.IOException;
import java.io.Reader;

public class PeekableReader extends Reader {
	private final Reader in;
	private boolean closed;
	private StringBuilder buffer;

	public PeekableReader(Reader in) {
		this.in = in;
		this.buffer = new StringBuilder();
		this.closed = false;
	}

	@Override
	public int read() throws IOException {
		int result;
		if (buffer.length() > 0) {
			result = buffer.charAt(0);
			buffer.delete(0, 1);
			return result;
		}
		return in.read();
	}

	@Override
	public int read(char[] cbuf, int off, int len) throws IOException {
		if (closed) {
			throw new IOException("Closed Reader!");
		}
		int read = 0;
		for (int i = 0; i < Math.min(len, buffer.length()); i++) {
			cbuf[off + i] = buffer.charAt(i);
			read++;
		}
		buffer.delete(0, read);
		if (read < len) {
			int num = in.read(cbuf, off + read, len - read);
			if (num == -1) {
				if (read == -1) {
					return -1;
				}
				return read;
			}
			read += num;
		}
		return read;
	}

	@Override
	public boolean ready() throws IOException {
		return buffer.length() > 0 || in.ready();
	}

	@Override
	public long skip(long n) throws IOException {
		long skipped = Math.min(n, buffer.length());
		buffer.delete(0, (int) skipped);
		if (skipped < n) {
			skipped += in.skip(n - skipped);
		}
		return skipped;
	}

	@Override
	public void close() throws IOException {
		in.close();
		this.closed = true;
	}
}
