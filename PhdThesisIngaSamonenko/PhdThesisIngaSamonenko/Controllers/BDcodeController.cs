using System.Web.Mvc;

namespace PhdThesisIngaSamonenko.Controllers
{
    public class BDcodeController : Controller
    {
        public ActionResult Index()
        {
            return File("BD.R", "text/plain");
        }

    }
}
