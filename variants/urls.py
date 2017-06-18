from django.conf.urls import url

from . import views

urlpatterns = [

    url(r'^variant_view$', views.variant_view, name='variant_view'),
    url(r'^Region_view$', views.Region_view, name='Region_view'),
    url(r'^Gene_view$', views.Gene_view, name='Gene_view'),
    url(r'^Transcript_view$', views.Transcript_view, name='Transcript_view'),
    url(r'^Home_view$', views.Home_view, name='Home_view'),
    url(r'^Advanced_Search$', views.Advanced_Search_view, name='Advanced_Search_view'),
    url(r'^Route$', views.Route, name='Route'),
    url(r'^ComplexQuery_view$', views.ComplexQuery_view, name='ComplexQuery_view'),
    url(r'^ComplexQuerydebug_view$', views.ComplexQuerydebug_view, name='ComplexQuerydebug_view')

]