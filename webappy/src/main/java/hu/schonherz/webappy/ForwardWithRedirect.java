package hu.schonherz.webappy;

import java.io.IOException;
import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * Servlet implementation class WithParam
 */
@WebServlet("/ForwardWithRedirect")
public class ForwardWithRedirect extends HttpServlet {
	private static final long serialVersionUID = 1L;
       
    /**
     * @see HttpServlet#HttpServlet()
     */
    public ForwardWithRedirect() {
        super();
        // TODO Auto-generated constructor stub
    }

	/**
	 * @see HttpServlet#doGet(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		String s = request.getParameter("parameter");
		String s2 = request.getParameter("parameter2");
		response.getWriter().write("doGet:" + s);
		System.out.println(s);
		System.out.println(s2);
		response.sendRedirect("WithParam?parameter="+s+"");
	}

	/**
	 * @see HttpServlet#doPost(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		String s = request.getParameter("parameter");
		String s2 = request.getParameter("parameter2");
		response.getWriter().write("doPost:" + s);
		System.out.println(s);
		System.out.println(s2);
		response.sendRedirect("WithParam");
		
	}

}
