clear all, close;

%�������� ������
k = 4;
m = [1 0 1 1];
l = -2: 1: 7;
d = 3;

g = [1 0 1 1];
e = [1 1 1 1 1 1 1];

% g = [1 1 0 1 1];
% e = [1 1 1 1 1 1 1 1];

%e = randi([0 1], 1, k + length(g) - 1);

%��������� ��������� �� ���� ������
[a,c] = code(g, m);
print("c(x)", c);
print("a(x)", a);

%������� ����� ��������� � ����� channel
print("e(x)",e);
b = channel(a, e);

%������� ����� ��������� �� ���� ��������
E = decode(g, b);
disp(E);

%�������������� �������
disp("for x^3+x+1 and k = " + k + " l = 3");
disp(dop(g, k, 3, d));

dmins_array = zeros(1, length(l));
for i = 1: length(l)
    dmins_array(i) = dmins(g, l(i) + k);
end

figure();
plot(l + k, dmins_array);
xlabel("����� �������������� �����");
ylabel("����������� ���������� ����");
grid on;

%---------------�������, ������������ ������ ��������� �������-------------

%������������� ������ ������
function [a,c] = code(g, m)
    r = length(g);
    factor = zeros(1, r);
    factor(1) = 1;

    multi = produce(m, factor);
    
    [~, c] = devide(multi,g);

    c = modul2(c);
    
    a = sum_(multi, c);
end

%������������� ������ ������
function b = channel(a, e)

    b = xor_(a, e);
    b = double(b);
    
    print("b(x)", b);
end

%������������� ������ ��������
function E = decode(g, b)

    [~, s] = devide(b, g);
    
    s = modul2(s);
    print("s(x)", s);

    E = comparison(s);
end

%����� ���� ������� ���� � ����� ������ (d-1)
function [result] = dop(g, k, l, d)

    %������������ ������ ���������
    message_length = k + l;
    [message_list,num]   = permitations(2, message_length);
    for i = 1: num
        message_list(i,:) = modul2(message_list(i,:));
    end
        
    %������������ ������ ������� ����
    codeword_list = zeros(num, message_len + degree(g));
    result = zeros(l + 1, degree(g) + message_length); 
    
    counter = 1;
    for i = 1: num
        word = code(g, message_list(i,:));
        codeword_list(i,:) = word;
        
        %����� ����������� ������� ����
        if (weigth(word) <= (d-1))
            result(counter,:) = word;
            counter = counter + 1;
        end
    end
end

%������� ���������� d_min ��� �������� ����� ���������
function result = dmins(poly, message_len)

    %������������ ��������� ���������
    [message_list,num]   = permitations(2, message_len);
    for i = 1: num
        message_list(i,:) = modul2(message_list(i,:));
    end

    %������������ ������� �����
    code_words = zeros(num, message_len + degree(poly));
    for i = 1: num
       [code_words(i,:), ~] = code(poly, message_list(i,:)); 
    end
    
    %������� ������������ ���� �������� ����� (!=0)
    code_word_weigth = zeros(1, num - 1);
    for i = 1: num - 1
        code_word_weigth(i) = weigth(code_words(i, :));
    end
    
    result = min(code_word_weigth);
end

%---------------------==������� ���������� � ����������==------------------

%������� ��������� ���� ���������
function result = produce(factor1, factor2)
    
    deg1 = length(factor1) - 1;
    deg2 = length(factor2) - 1;
    deg_result = deg1 + deg2;
    
    result = zeros(1, deg_result + 1);

    for i = 1: length(factor1)
       for j = 1: length(factor2)
           result(i + j - 1) = result(i + j - 1) + ( factor1(i) * factor2(j) );
       end
    end
    
    for i = 1: length(result)
        result(i) = mod2(result(i)); 
    end
end

%������� ������� ���� ��������� � ��������
function [result, remainder] = devide(first, second)

    first  = polyCut(first);
    second = polyCut(second);

    deg1 = degree(first);
    deg2 = degree(second);
    deg  = deg1 - deg2;

    remainder = zeros(1, deg1 + 1);
    
    if (deg1 < deg2)
        %������� ������ ��������
        result    = zeros(1, deg  + 1);
        remainder = first;
    else
        
        tmp_deg = deg1 - deg2;
        result  = zeros(1, tmp_deg  + 1);
        result(length(result) - tmp_deg) = 1;
        
        while (tmp_deg >= 0)

            ff = produce(second, result);
            remainder = first - ff;
            remainder = modul2(remainder);
            tmp_deg = degree(remainder) - deg2;
            if (degree(remainder) >= deg2)
                result(length(result) - tmp_deg) = 1;
            end
        end
    end
end

%������� ������������� xor-a ���� ���������
function result = xor_(poly1, poly2)
    result = zeros(1, length(poly1));
    
    for i = 1: length(result)
        if(poly2(i) == 1)
            result(i) = ~(poly1(i));
        else
            result(i) = poly1(i);
        end
    end
end

%������� ���������� ������� ��������
function result = degree(poly)
    i = 1;
    while ((poly(i) ~= 1) && ( i < length(poly)))
        i = i + 1;
    end
    result = length(poly) - i;
end

%������� ���������� �������� � ����������� ����
function result = polyCut(poly)
    deg = degree(poly);
    result = zeros(1, deg + 1);
    
    counter = 1;
    for i = length(poly) - deg: length(poly)
        result(counter) = poly(i);
        counter = counter + 1;
    end
end

%������� �������� ���� ��������� (�������� ������ ������)
function result = sum_(poly1, poly2)

    len1 = length(poly1);
    len2 = length(poly2);
    d_len = len1 - len2;
    
    if (d_len == 0)
        result = poly1 + poly2;
    else
        if (d_len < 0)
            tmp = zeros(1, len1 + abs(d_len));
            counter = 1;
            for i = abs(d_len): length(tmp)
                tmp(i) = poly1(counter);
                counter = counter + 1;
            end
            result = tmp + poly2;
        else
            tmp = zeros(1, len2 + abs(d_len));
            counter = 1;
            for i = (abs(d_len) + 1): length(tmp)
                tmp(i) = poly2(counter);
                counter = counter + 1;
            end
            result = poly1 + tmp;
        end
    end
end

%--------------------------��������������� �������-------------------------

%������� ������ ��������
function print(string, array)
    disp(string + " = ");
    disp(array);
end

%������� ���������� ������������� �������� �� ������
function array = modul2(vector)
    array = zeros(1, length(vector));
    for i = 1: length(vector)
        array(i) = mod2(vector(i));
    end
end

%������� ��������� ������� � 0
function comp = comparison(array)
    comp = 0;
    for i = 1: length(array)
       comp = comp + array(i); 
    end
    
    if(comp ~= 0)
        comp = 1;
    else
        comp = 0;
    end  
end

%������� ���������� ���� �������
function w = weigth(vector)
    w = 0;
    for i = 1: length(vector)
        w = w + vector(i);
    end
end

%������� ���������� ����� �� ������
function r = mod2(num)
    t = int8(num / 2);
    r = num - (2 * t);
    if(r < 0)
        r = r + 2;
    end
end

%������� ���������� ���������� �� n �� k
function [result, num] = permitations(n, k)

    num = n^k;
    result = zeros(num, k);
    
    for i = 1:  num
        result(i,:) = getSequence(i, n, k);
    end
end

function result = getSequence(sequence, n, k)
    result = zeros(1, k);
    while (k > 0)
        result(k) = rem(sequence,n);
        sequence = int8(sequence / n);
        k = k - 1;
    end
end