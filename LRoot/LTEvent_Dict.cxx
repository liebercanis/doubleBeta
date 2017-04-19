// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LTEvent_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "LTEvent.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LTEvent(void *p = 0);
   static void *newArray_LTEvent(Long_t size, void *p);
   static void delete_LTEvent(void *p);
   static void deleteArray_LTEvent(void *p);
   static void destruct_LTEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LTEvent*)
   {
      ::LTEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LTEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LTEvent", ::LTEvent::Class_Version(), "LTEvent.hxx", 22,
                  typeid(::LTEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LTEvent::Dictionary, isa_proxy, 4,
                  sizeof(::LTEvent) );
      instance.SetNew(&new_LTEvent);
      instance.SetNewArray(&newArray_LTEvent);
      instance.SetDelete(&delete_LTEvent);
      instance.SetDeleteArray(&deleteArray_LTEvent);
      instance.SetDestructor(&destruct_LTEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LTEvent*)
   {
      return GenerateInitInstanceLocal((::LTEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LTEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LTEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LTEvent::Class_Name()
{
   return "LTEvent";
}

//______________________________________________________________________________
const char *LTEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LTEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LTEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LTEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LTEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class LTEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LTEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(LTEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LTEvent(void *p) {
      return  p ? new(p) ::LTEvent : new ::LTEvent;
   }
   static void *newArray_LTEvent(Long_t nElements, void *p) {
      return p ? new(p) ::LTEvent[nElements] : new ::LTEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_LTEvent(void *p) {
      delete ((::LTEvent*)p);
   }
   static void deleteArray_LTEvent(void *p) {
      delete [] ((::LTEvent*)p);
   }
   static void destruct_LTEvent(void *p) {
      typedef ::LTEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LTEvent

namespace ROOT {
   static TClass *vectorlELTTrajectgR_Dictionary();
   static void vectorlELTTrajectgR_TClassManip(TClass*);
   static void *new_vectorlELTTrajectgR(void *p = 0);
   static void *newArray_vectorlELTTrajectgR(Long_t size, void *p);
   static void delete_vectorlELTTrajectgR(void *p);
   static void deleteArray_vectorlELTTrajectgR(void *p);
   static void destruct_vectorlELTTrajectgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LTTraject>*)
   {
      vector<LTTraject> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LTTraject>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LTTraject>", -2, "vector", 210,
                  typeid(vector<LTTraject>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELTTrajectgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<LTTraject>) );
      instance.SetNew(&new_vectorlELTTrajectgR);
      instance.SetNewArray(&newArray_vectorlELTTrajectgR);
      instance.SetDelete(&delete_vectorlELTTrajectgR);
      instance.SetDeleteArray(&deleteArray_vectorlELTTrajectgR);
      instance.SetDestructor(&destruct_vectorlELTTrajectgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LTTraject> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<LTTraject>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELTTrajectgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LTTraject>*)0x0)->GetClass();
      vectorlELTTrajectgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELTTrajectgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELTTrajectgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LTTraject> : new vector<LTTraject>;
   }
   static void *newArray_vectorlELTTrajectgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LTTraject>[nElements] : new vector<LTTraject>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELTTrajectgR(void *p) {
      delete ((vector<LTTraject>*)p);
   }
   static void deleteArray_vectorlELTTrajectgR(void *p) {
      delete [] ((vector<LTTraject>*)p);
   }
   static void destruct_vectorlELTTrajectgR(void *p) {
      typedef vector<LTTraject> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LTTraject>

namespace ROOT {
   static TClass *vectorlELTPVertexgR_Dictionary();
   static void vectorlELTPVertexgR_TClassManip(TClass*);
   static void *new_vectorlELTPVertexgR(void *p = 0);
   static void *newArray_vectorlELTPVertexgR(Long_t size, void *p);
   static void delete_vectorlELTPVertexgR(void *p);
   static void deleteArray_vectorlELTPVertexgR(void *p);
   static void destruct_vectorlELTPVertexgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LTPVertex>*)
   {
      vector<LTPVertex> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LTPVertex>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LTPVertex>", -2, "vector", 210,
                  typeid(vector<LTPVertex>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELTPVertexgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<LTPVertex>) );
      instance.SetNew(&new_vectorlELTPVertexgR);
      instance.SetNewArray(&newArray_vectorlELTPVertexgR);
      instance.SetDelete(&delete_vectorlELTPVertexgR);
      instance.SetDeleteArray(&deleteArray_vectorlELTPVertexgR);
      instance.SetDestructor(&destruct_vectorlELTPVertexgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LTPVertex> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<LTPVertex>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELTPVertexgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LTPVertex>*)0x0)->GetClass();
      vectorlELTPVertexgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELTPVertexgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELTPVertexgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LTPVertex> : new vector<LTPVertex>;
   }
   static void *newArray_vectorlELTPVertexgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LTPVertex>[nElements] : new vector<LTPVertex>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELTPVertexgR(void *p) {
      delete ((vector<LTPVertex>*)p);
   }
   static void deleteArray_vectorlELTPVertexgR(void *p) {
      delete [] ((vector<LTPVertex>*)p);
   }
   static void destruct_vectorlELTPVertexgR(void *p) {
      typedef vector<LTPVertex> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LTPVertex>

namespace {
  void TriggerDictionaryInitialization_LTEvent_Dict_Impl() {
    static const char* headers[] = {
"LTEvent.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/usr/local/root-6.08.00-build/include",
"/home/gold/doubleBeta/LRoot/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LTEvent_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$LTEvent.hxx")))  LTPVertex;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$LTEvent.hxx")))  LTEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LTEvent_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "LTEvent.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"LTEvent", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LTEvent_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LTEvent_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LTEvent_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LTEvent_Dict() {
  TriggerDictionaryInitialization_LTEvent_Dict_Impl();
}
