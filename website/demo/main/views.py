from django.shortcuts import render 
from django.http import HttpResponse
from matplotlib.pyplot import *
import matplotlib
matplotlib.use('Agg') 
from matplotlib.backends.backend_agg import FigureCanvasAgg
from numpy import *
import base64
from io import BytesIO
from .models import Tutorial
# Create your views here.


def homepage(request):
    # fig, ax = subplots()
    # ax.plot(range(10))
    # response = HttpResponse(content_type = 'image/png')
    # canvas = FigureCanvasAgg(fig)
    # canvas.print_png(response)
    # return response
    # return HttpResponse("Wow this is an awesome toturial")
    # if request.method=="POST":
    #     form = request.POST.dict()
    #     print(form)
    #     return render(request=request, 
    #         template_name="main/home.html",
    #         context = {"toturial":Tutorial.objects.all})

    return render(request=request, 
            template_name="main/home.html")


def SNRT(request):
    if request.method=="POST":
        form = request.POST.dict()
        print(form)
        return render(request=request, 
            template_name="main/SNRT.html",
            context = {"img":t_snr(form)})

    return render(request=request, 
            template_name="main/SNRT.html",
            context = {"toturial":Tutorial.objects.all})

def SNRNs(request):
    if request.method=="POST":
        form = request.POST.dict()
        print(form)
        return render(request=request, 
            template_name="main/SNRNs.html",
            context = {"img":step_snr(form)})

    return render(request=request, 
            template_name="main/SNRNs.html",
            context = {"toturial":Tutorial.objects.all})

def SNRNa(request):
    if request.method=="POST":
        form = request.POST.dict()
        print(form)

        return render(request=request, 
            template_name="main/SNRNa.html",
            context = {"img":na_snr(form)})
    return render(request=request, 
            template_name="main/SNRNa.html",
            context = {"toturial":Tutorial.objects.all})


def g_alpha(request):
    if request.method=="POST":
        form = request.POST.dict()
        return render(request=request,
            template_name="main/g_alpha.html",
            context = {"limit":g_limit(form)})
    return render(request=request, 
            template_name="main/g_alpha.html",
            context = {"toturial":Tutorial.objects.all}) 


g_KSVZ = -0.94
alpha = 1/137

def t_snr(di):
    scan_range = eval(di["scan_range"])
    step_size = eval(di["step_size"])
    delta_v = eval(di["delta_v"])
    rho_a = eval(di["rho_a"])
    f = eval(di["f"])
    B = eval(di["B"])
    V = eval(di["V"])
    T = eval(di["T"])
    Q = eval(di["Q"])
    c = eval(di["c"])
    s = eval(di["s"])
    Ns = scan_range / step_size
    m = 4.15e-6 * f
    g_target_13 = 3e-13 * pi *0.62 * 1e7/ (alpha*m)
    g_target_14  = 9.979108e-15 * pi *0.62 * 1e7/ (alpha*m)

    SNR = linspace(0,10,1000)
    t = []
    up  = (g_KSVZ/0.97)**2*(rho_a/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B**2)*(V)* sqrt(Ns)
    low = SNR * sqrt(delta_v)
    t= (0.0017 * up / low)**-2
    try:
        close()
        title("time to reach KSVZ")
        xlabel("SNR")
        ylabel("t")
        plot(SNR,t)
        xlim(0,10)
        grid()
        fill_between(SNR,t,0)
    except:
        print("error")
        return None
    buffer = BytesIO()
    savefig(buffer)
    close()
    plot_data = buffer.getvalue()
    imb = base64.b64encode(plot_data)
    ims = imb.decode()
    imd = "data:image/png;base64,"+ims
    return imd 

def step_snr(di):
    scan_range = eval(di["scan_range"])
    Na = eval(di["Na"])
    delta_v = eval(di["delta_v"])
    rho_a = eval(di["rho_a"])
    f = eval(di["f"])
    B = eval(di["B"])
    V = eval(di["V"])
    T = eval(di["T"])
    Q = eval(di["Q"])
    c = eval(di["c"])
    s = eval(di["s"])
    m = 4.15e-6 * f
    g_target_13 = 3e-13 * pi *0.62 * 1e7/ (alpha*m)
    g_target_14  = 9.979108e-15 * pi *0.62 * 1e7/ (alpha*m)

    SNR = linspace(-0.1,10,1000)
    step = []
    up  = (g_KSVZ/0.97)**2*(rho_a/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B**2)*(V)* sqrt(scan_range * Na)
    low = SNR * delta_v
    KSVZ_step = (0.0017 * up / low)**2
    up  = (g_target_13/0.97)**2*(rho_a/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B**2)*(V)* sqrt(scan_range * Na)
    step_13 = (0.0017 * up / low)**2
    up  = (g_target_14/0.97)**2*(rho_a/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B**2)*(V)* sqrt(scan_range * Na)
    step_14 = (0.0017 * up / low)**2
    try:
        close()
        title("SNR vs step")
        xlabel("SNR")
        ylabel("step")
        plot(SNR,KSVZ_step,label="KSVZ", alpha = 0.8,color="#2b79ff")
        fill_between(SNR,KSVZ_step,0, alpha = 0.8,color="#2b79ff")
        plot(SNR,step_13,label="$g_{\\alpha\\gamma\\gamma}$ = 3e-13", alpha =  0.8,color="#54ffaa")

        plot(SNR,step_14,label="$g_{\\alpha\\gamma\\gamma}$ = 9.97e-15", alpha =  0.8,color="#ff5454")
        fill_between(SNR,step_14,KSVZ_step, alpha =  0.8,color="#ff5454")
        fill_between(SNR,step_13,step_14, alpha =  0.8,color="#54ffaa")
        legend()
        xlim(0,10)
        grid()
        yscale("log")
    except:
        return None
    buffer = BytesIO()
    savefig(buffer)
    close()
    plot_data = buffer.getvalue()
    imb = base64.b64encode(plot_data)
    ims = imb.decode()
    imd = "data:image/png;base64,"+ims
    return imd 

def na_snr(di):
    scan_range = eval(di["scan_range"])
    step_size = eval(di["step_size"])
    delta_v = eval(di["delta_v"])
    rho_a = eval(di["rho_a"])
    f = eval(di["f"])
    B = eval(di["B"])
    V = eval(di["V"])
    T = eval(di["T"])
    Q = eval(di["Q"])
    c = eval(di["c"])
    s = eval(di["s"])
    Ns = scan_range / step_size
    m = 4.15e-6 * f
    g_target_13 = 3e-13 * pi *0.62 * 1e7/ (alpha*m)
    g_target_14  = 9.979108e-15 * pi *0.62 * 1e7/ (alpha*m)

    SNR = linspace(0,10,1000)
    na = []
    up  = (g_KSVZ/0.97)**2*(rho_a/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B**2)*(V)* sqrt(Ns)
    low = SNR * delta_v
    t_KSVZ = (0.0017 * up / low)**-2

    up  = (g_target_13/0.97)**2*(rho_a/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B**2)*(V)* sqrt(Ns)
    t_13 = (0.0017 * up / low)**-2

    up  = (g_target_14/0.97)**2*(rho_a/0.45)*(f/1)*(Q/10_000)*(c/0.5)*(1-2*s)/(1-s)*(B**2)*(V)* sqrt(Ns)
    t_14 = (0.0017 * up / low)**-2
    try:
        close()
        title("SNR vs Na")
        xlabel("SNR")
        ylabel("Na")
        plot(SNR,t_KSVZ, alpha = 0.8,color="#2b79ff",label="KSVZ")
        plot(SNR,t_13, alpha =  0.8,color="#54ffaa",label="$g_{\\alpha\\gamma\\gamma}$ = 3e-13")
        plot(SNR,t_14, alpha =  0.8,color="#ff5454",label="$g_{\\alpha\\gamma\\gamma}$ = 9.97e-15")

        fill_between(SNR,t_13,0, alpha =  0.8,color="#54ffaa")
        fill_between(SNR,t_14,t_13, alpha =  0.8,color="#ff5454")
        fill_between(SNR,t_KSVZ,t_14, alpha =  0.8,color="#2b79ff")
        legend()
        xlim(0,10)
        grid()
        yscale("log")
    except:
        return None
    buffer = BytesIO()
    savefig(buffer)
    close()
    plot_data = buffer.getvalue()
    imb = base64.b64encode(plot_data)
    ims = imb.decode()
    imd = "data:image/png;base64,"+ims
    return imd 

def g_limit(di):
    t = eval(di["time"])
    delta_v = eval(di["delta_v"])
    rho_a = eval(di["rho_a"]) * 1e15
    f = eval(di["f"]) * 1e9
    B = eval(di["B"]) 
    V = eval(di["V"]) * 1e-3
    T = eval(di["T"])
    Q = eval(di["Q"])
    c = eval(di["c"])
    s = eval(di["s"])
    sigma = eval(di["sigma"])
    beta = eval(di["beta"])
    k = 1.38e-23
    h  = 6.626e-34 # J/k
    hbar  = 4.135e-15/(2*pi)
    mu_0 = 4*pi  * 1e-7
    w_x = 2*pi * f
    m_a = hbar * w_x
    C  = 3e8
    sigma_1 = h * f * sqrt(delta_v/t)
    sigma_2 = k * T * sqrt(delta_v/t)
    shift = (hbar * C )**3 * rho_a / (m_a**2) * c * (beta)/(beta+1) * Q * w_x * V * (1/mu_0)* B**2
    first_line = "{0:e} : Thermal noise (kb T)".format(sqrt(sigma * sigma_2 / shift)*1e9 )
    second_line = "{0:e} : Quantum noise (h f)".format(sqrt(sigma * sigma_1 / shift)*1e9 )
    return [first_line,second_line]