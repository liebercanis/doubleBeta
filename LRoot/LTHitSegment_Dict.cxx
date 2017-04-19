// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LTHitSegment_Dict

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
#include "LTHitSegment.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LTHitSegment(void *p = 0);
   static void *newArray_LTHitSegment(Long_t size, void *p);
   static void delete_LTHitSegment(void *p);
   static void deleteArray_LTHitSegment(void *p);
   static void destruct_LTHitSegment(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LTHitSegment*)
   {
      ::LTHitSegment *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LTHitSegment >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LTHitSegment", ::LTHitSegment::Class_Version(), "LTHitSegment.hxx", 26,
                  typeid(::LTHitSegment), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LTHitSegment::Dictionary, isa_proxy, 4,
                  sizeof(::LTHitSegment) );
      instance.SetNew(&new_LTHitSegment);
      instance.SetNewArray(&newArray_LTHitSegment);
      instance.SetDelete(&delete_LTHitSegment);
      instance.SetDeleteArray(&deleteArray_LTHitSegment);
      instance.SetDestructor(&destruct_LTHitSegment);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LTHitSegment*)
   {
      return GenerateInitInstanceLocal((::LTHitSegment*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LTHitSegment*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LTHitSegment::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LTHitSegment::Class_Name()
{
   return "LTHitSegment";
}

//______________________________________________________________________________
const char *LTHitSegment::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTHitSegment*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LTHitSegment::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTHitSegment*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LTHitSegment::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTHitSegment*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LTHitSegment::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTHitSegment*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LTHitSegment::Streamer(TBuffer &R__b)
{
   // Stream an object of class LTHitSegment.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LTHitSegment::Class(),this);
   } else {
      R__b.WriteClassBuffer(LTHitSegment::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LTHitSegment(void *p) {
      return  p ? new(p) ::LTHitSegment : new ::LTHitSegment;
   }
   static void *newArray_LTHitSegment(Long_t nElements, void *p) {
      return p ? new(p) ::LTHitSegment[nElements] : new ::LTHitSegment[nElements];
   }
   // Wrapper around operator delete
   static void delete_LTHitSegment(void *p) {
      delete ((::LTHitSegment*)p);
   }
   static void deleteArray_LTHitSegment(void *p) {
      delete [] ((::LTHitSegment*)p);
   }
   static void destruct_LTHitSegment(void *p) {
      typedef ::LTHitSegment current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LTHitSegment

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace {
  void TriggerDictionaryInitialization_LTHitSegment_Dict_Impl() {
    static const char* headers[] = {
"LTHitSegment.hxx",
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
#line 1 "LTHitSegment_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$LTHitSegment.hxx")))  LTHitSegment;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LTHitSegment_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "LTHitSegment.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"LTHitSegment", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LTHitSegment_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LTHitSegment_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LTHitSegment_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LTHitSegment_Dict() {
  TriggerDictionaryInitialization_LTHitSegment_Dict_Impl();
}
