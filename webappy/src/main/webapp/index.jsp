<html>
<head><title>Hello World</title></head>
<body>
	<form action="ForwardWithParam" method="get">
		<input type="text" name="parameter">
		<input type="text" name="parameter2">
		<input type="submit">
	</form>
	<% String s = (String) request.getAttribute("result");	%>
	
</body>
</html>